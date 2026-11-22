*! version 1.0.5  21nov2026 nwlsm
* Copyright (C) 2026 Thomas U. Grund
* This file is part of nwlsm, a Stata package for latent space models for networks. It follows the latent space model of
* Hoff, Raftery and Handcock (2002) as implemented in the R package latentnet (Krivitsky and Handcock) and is written independently.
* nwlsm is free software: you can redistribute it and/or modify it under the terms of the GNU General
* Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
* that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.

capture program drop nwlsm
program nwlsm, eclass
	version 14
	gettoken sub rest : 0, parse(" ,")
	if "`sub'" == "version" {
		nwlsm_version
		exit
	}
	if "`sub'" == "" | substr(trim(`"`0'"'), 1, 1) == "," {
		if "`e(cmd)'" != "nwlsm" error 301
		_nwlsm_display
		exit
	}
	syntax anything(name=netname) [, Dim(integer 2) EDGEcov(string) ABSdiff(string) NSample(integer 4000) BURNin(integer 10000) ///
		INTerval(integer 10) SEED(string) PILOT(integer 4) PLEN(integer 500) RANDOMstart SAVEDraws(string) ///
		SENDER RECEIVER SOCiality SENDERVar(real 1) SENDERDf(real 3) RECEIVERVar(real 1) RECEIVERDf(real 3) SOCIALITYVar(real 1) SOCIALITYDf(real 3) BILinear MKLStarts(integer -1) CHAins(integer 1) THReads(integer 0) CLUsters(integer 0) FAMily(string) TRIals(string) DISPVar(real 1) DISPDf(real 3) ///
		INFluence(string) NEIghbors(integer 0) WEIght(string) OBVar(real 100) OSIGma(numlist min=2 max=2 >0)]
	if `dim' < 1 | `dim' > 10 {
		di as error "dim() must be an integer between 1 and 10"
		error 198
	}
	local __fam 0
	local __famname "bernoulli"
	if `"`family'"' != "" {
		local __fl = lower(`"`family'"')
		if inlist("`__fl'", "bernoulli", "binary") {
			local __fam 0
		}
		else if "`__fl'" == "poisson" {
			local __fam 1
			local __famname "poisson"
		}
		else if "`__fl'" == "binomial" {
			local __fam 2
			local __famname "binomial"
		}
		else if inlist("`__fl'", "normal", "gaussian") {
			local __fam 3
			local __famname "normal"
		}
		else {
			di as error "{p 0 0 2}family() must be bernoulli (the default), poisson, binomial or normal.{p_end}"
			error 198
		}
	}
	if `"`trials'"' != "" & `__fam' != 2 {
		di as error "trials() belongs to family(binomial)"
		error 198
	}
	if `__fam' == 2 & `"`trials'"' == "" {
		di as error "{p 0 0 2}family(binomial) needs trials(): a number, a network or a matrix with the numbers of trials of the dyads.{p_end}"
		error 198
	}
	if `__fam' == 3 & (`dispvar' <= 0 | `dispdf' <= 0) {
		di as error "dispvar() and dispdf() (the prior scale and degrees of freedom of the dispersion) must be positive"
		error 198
	}
	if `nsample' < 10 | `burnin' < 0 | `interval' < 1 | `pilot' < 0 | `plen' < 10 {
		di as error "{p 0 0 2}nsample() must be at least 10, burnin() not negative, interval() at least 1, pilot() not negative and plen() at least 10.{p_end}"
		error 198
	}
	if `sendervar' <= 0 | `receivervar' <= 0 | `socialityvar' <= 0 | `senderdf' <= 0 | `receiverdf' <= 0 | `socialitydf' <= 0 {
		di as error "the prior scales and degrees of freedom of the random effects must be positive"
		error 198
	}
	if `mklstarts' != -1 & `mklstarts' < 1 {
		di as error "mklstarts() must be at least 1"
		error 198
	}
	* several MKL starts by default when random effects are in the model (their MKL estimates are weakly identified and one start can end in a poor local optimum)
	if `mklstarts' == -1 local mklstarts = cond("`sender'`receiver'`sociality'" != "", 4, 1)
	if `chains' < 1 | `chains' > 64 {
		di as error "chains() must be an integer between 1 and 64"
		error 198
	}
	if `threads' < 0 | `threads' > 256 {
		di as error "threads() must be between 1 and 256 (default: 4 threads for one chain, 2 per chain for several chains, at most one per core of the computer)"
		error 198
	}
	// Default number of threads (threads(0)): one chain gets 4 threads for its per-dyad passes, several chains 2 threads per chain, at most the processors of the machine.
	// More helper threads per chain make the per-dyad passes slower, not faster (measured at 200 actors); the estimates do not depend on the number of threads.
	local __thr = `threads'
	if `__thr' == 0 {
		local __thr = min(c(processors_mach), cond(`chains' == 1, 4, 2 * `chains'))
		if `__thr' < 1 local __thr 1
	}
	if `clusters' != 0 & (`clusters' < 2 | `clusters' > 6) {
		di as error "clusters() must be an integer between 2 and 6"
		error 198
	}
	if ceil(`nsample' / `chains') < 10 {
		di as error "nsample() must give at least 10 draws for each of the chains()"
		error 198
	}
	if `"`seed'"' != "" {
		capture confirm integer number `seed'
		if _rc | real(`"`seed'"') < 0 | real(`"`seed'"') > 2147483647 {
			di as error "seed() must be an integer between 0 and 2147483647"
			error 198
		}
	}
	local __io 0
	if `"`influence'"' != "" {
		local __io 1
		if "`bilinear'" != "" {
			di as error "influence() needs distances between the latent positions: it cannot be combined with bilinear"
			error 198
		}
		if `neighbors' < 0 {
			di as error "neighbors() must be a positive integer (the default, 0, uses all other actors)"
			error 198
		}
		if `obvar' <= 0 {
			di as error "obvar() must be positive"
			error 198
		}
		local __wt 0
		local __wtname "invdist"
		if `"`weight'"' != "" {
			local __wl = lower(`"`weight'"')
			if inlist("`__wl'", "invdist", "inverse") local __wt 0
			else if inlist("`__wl'", "expdist", "exp") {
				local __wt 1
				local __wtname "expdist"
			}
			else {
				di as error "{p 0 0 2}weight() must be invdist (weights proportional to 1/distance, as in Sweet and Adhikari 2020) or expdist (weights proportional to exp(-distance)).{p_end}"
				error 198
			}
		}
		local __osa 0.01
		local __osb 0.01
		if "`osigma'" != "" {
			local __osa : word 1 of `osigma'
			local __osb : word 2 of `osigma'
		}
	}
	else if `neighbors' != 0 | `"`weight'"' != "" | "`osigma'" != "" | `obvar' != 100 {
		di as error "neighbors(), weight(), obvar() and osigma() belong to influence()"
		error 198
	}
	_nwlsm_libcheck
	set more off
	_nwsyntax `netname', max(1)
	if "`is2mode'" == "true" {
		di as error "nwlsm fits one-mode networks; `netname' is a two-mode network."
		error 198
	}
	if "`valued'" == "true" & `__fam' == 0 {
		mata: st_local("__valok", strofreal(`netobj'->get_arcs_count() == 0 | (`netobj'->get_maximum() <= 1 & `netobj'->get_minimum() >= 0 & `netobj'->get_arcs_sum() == `netobj'->get_arcs_count())))
		if !`__valok' {
			di as error "{p 0 0 2}`netname' is a valued network: use family(poisson), family(binomial) with trials(), or family(normal) to model the values, or dichotomise it first with {bf:nwreplace}.{p_end}"
			error 198
		}
	}
	if `nodes' < 3 {
		di as error "the network needs at least 3 actors"
		error 198
	}
	local __dir = ("`directed'" == "true")
	local __onames
	if `__io' {
		capture unab __ivars : `influence'
		if _rc {
			di as error "influence(): not a list of variables in the current dataset"
			error 111
		}
		local __ni : word count `__ivars'
		if `__ni' < 2 {
			di as error "{p 0 0 2}influence() needs the outcome variable and the lagged outcome variable (the outcome measured before the network), then any other regressors.{p_end}"
			error 198
		}
		local __ivd : list dups __ivars
		if "`__ivd'" != "" {
			di as error "influence(): the same variable is listed twice"
			error 198
		}
		confirm numeric variable `__ivars'
		if _N != `nodes' {
			di as error "{p 0 0 2}influence() needs one observation per actor in the current dataset, in the order of the actors of {bf:`netname'} (the dataset has `=_N' observations, the network `nodes' actors).{p_end}"
			error 198
		}
		if `__ni' > 39 {
			di as error "influence(): at most 37 regressors besides the lagged outcome"
			error 198
		}
		foreach __v of local __ivars {
			mata: st_local("__hm", strofreal(hasmissing(st_data(., "`__v'"))))
			if `__hm' {
				di as error "influence(): variable {bf:`__v'} has missing values"
				error 416
			}
		}
		local __ov : word 1 of `__ivars'
		local __lv : word 2 of `__ivars'
		local __onames "_cons `=subinstr("`__ivars'", "`__ov' ", "", 1)' influence"
		capture mata: mata drop __nwlsm_in
		mata: __nwlsm_in = st_data(., ("`__ov'", "`__lv'")), J(`nodes', 1, 1), st_data(., tokens("`=subinstr("`__ivars'", "`__ov' ", "", 1)'"))
	}
	if !`__dir' & ("`sender'" != "" | "`receiver'" != "") {
		di as error "{p 0 0 2}sender and receiver effects need a directed network; use {bf:sociality} for the undirected network {bf:`netname'}.{p_end}"
		error 198
	}
	if "`sociality'" != "" & ("`sender'" != "" | "`receiver'" != "") & !`__dir' {
		di as error "sociality cannot be combined with sender or receiver effects in an undirected network"
		error 198
	}
	local __opt "`=("`sender'" != "")' `=("`receiver'" != "")' `=("`sociality'" != "")' `=("`bilinear'" != "")' `sendervar' `senderdf' `receivervar' `receiverdf' `socialityvar' `socialitydf' `mklstarts' `chains' `__thr' `clusters' `__fam' `dispvar' `dispdf'`=cond(`__io', " `neighbors' `__wt' `obvar' `__osa' `__osb'", "")'"
	capture mata: mata drop __nwlsm_y
	capture mata: mata drop __nwlsm_X
	capture mata: mata drop __nwlsm_tr
	if `__fam' == 0 mata: __nwlsm_y = lsm_ymat(*(`netobj'->get_matrix()))
	else mata: __nwlsm_y = *(`netobj'->get_matrix())
	mata: __nwlsm_tr = J(0, 0, .)
	mata: __nwlsm_X = J(0, `nodes', .)
	local __names "_cons"
	foreach __c of local edgecov {
		capture _nwsyntax `__c', max(1) other(__cv)
		if !_rc {
			if "`__cvis2mode'" == "true" | `__cvnodes' != `nodes' {
				capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
				di as error "{p 0 0 2}the edge covariate {bf:`__c'} must be a one-mode network of the same `nodes' actors as {bf:`netname'}.{p_end}"
				error 198
			}
			mata: __nwlsm_X = __nwlsm_X \ *(`__cvnetobj'->get_matrix())
		}
		else {
			capture confirm matrix `__c'
			if _rc {
				capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
				di as error "{p 0 0 2}edgecov(): {bf:`__c'} is neither a network nor a Stata matrix.{p_end}"
				error 198
			}
			if rowsof(`__c') != `nodes' | colsof(`__c') != `nodes' {
				capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
				di as error "{p 0 0 2}the matrix {bf:`__c'} must have `nodes' rows and columns, one for each actor of {bf:`netname'}.{p_end}"
				error 198
			}
			mata: __nwlsm_X = __nwlsm_X \ st_matrix("`__c'")
		}
		local __names `__names' edgecov_`__c'
	}
	if "`absdiff'" != "" {
		capture unab absdiff : `absdiff'
		if _rc {
			capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
			di as error "absdiff(): not a list of variables in the current dataset"
			error 111
		}
		capture confirm numeric variable `absdiff'
		if _rc {
			capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
			di as error "absdiff(): the variables must be numeric"
			error 109
		}
		if _N != `nodes' {
			capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
			di as error "{p 0 0 2}absdiff() needs one observation per actor in the current dataset, in the order of the actors of {bf:`netname'} (the dataset has `=_N' observations, the network `nodes' actors).{p_end}"
			error 198
		}
		foreach __v of local absdiff {
			mata: st_local("__hm", strofreal(hasmissing(st_data(., "`__v'"))))
			if `__hm' {
				capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
				di as error "absdiff(): variable {bf:`__v'} has missing values"
				error 416
			}
			mata: __nwlsm_X = __nwlsm_X \ lsm_node_covs("`__v'", `nodes')
			local __names `__names' absdiff_`__v'
		}
	}
	* the diagonal of a one-mode matrix is not a dyad: a missing diagonal (as in a network matrix) is set to 0, missing off-diagonal values stay an error
	mata: st_local("__hm", strofreal(rows(__nwlsm_X) > 0 ? any((__nwlsm_X :>= .) :* (1 :- (J(rows(__nwlsm_X) / `nodes', 1, 1) # I(`nodes')))) : 0))
	if !`__hm' mata: __nwlsm_X = editmissing(__nwlsm_X, 0)
	if `__hm' {
		capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
		di as error "edge covariates must not have missing values"
		error 416
	}
	local __nc : word count `__names'
	local __dp : list dups __names
	if "`__dp'" != "" {
		capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
		di as error "the same edge covariate is listed twice"
		error 198
	}
	if `__fam' > 0 {
		if `__fam' == 2 {
			capture confirm number `trials'
			if !_rc {
				if `trials' < 1 | `trials' != floor(`trials') {
					capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
					di as error "trials() must be a positive integer, a network or a matrix of positive integers"
					error 198
				}
				mata: __nwlsm_tr = J(`nodes', `nodes', `trials')
			}
			else {
				capture _nwsyntax `trials', max(1) other(__tr)
				if !_rc {
					if "`__tris2mode'" == "true" | `__trnodes' != `nodes' {
						capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
						di as error "{p 0 0 2}trials(): the network {bf:`trials'} must be a one-mode network of the same `nodes' actors as {bf:`netname'}.{p_end}"
						error 198
					}
					mata: __nwlsm_tr = *(`__trnetobj'->get_matrix())
				}
				else {
					capture confirm matrix `trials'
					if _rc {
						capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
						di as error "{p 0 0 2}trials(): {bf:`trials'} is neither a number, a network nor a Stata matrix.{p_end}"
						error 198
					}
					if rowsof(`trials') != `nodes' | colsof(`trials') != `nodes' {
						capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
						di as error "{p 0 0 2}the matrix {bf:`trials'} must have `nodes' rows and columns, one for each actor of {bf:`netname'}.{p_end}"
						error 198
					}
					mata: __nwlsm_tr = st_matrix("`trials'")
				}
			}
		}
		mata: st_local("__vc", strofreal(lsm_checkvalues(__nwlsm_y, `__dir', `__fam', __nwlsm_tr)))
		if `__vc' > 0 {
			capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr
			if `__vc' == 1 di as error "{p 0 0 2}family(`__famname') needs non-negative integer values on the dyads of {bf:`netname'}.{p_end}"
			if `__vc' == 2 di as error "{p 0 0 2}a value of {bf:`netname'} is larger than the number of trials of its dyad.{p_end}"
			if `__vc' == 4 di as error "{p 0 0 2}trials() must be a positive integer for every observed dyad.{p_end}"
			error 198
		}
	}
	* degenerate network: no observed tie at all or every observed dyad a tie (binary); all observed values 0 (poisson, binomial); the diagonal is ignored
	mata: st_local("__degen", strofreal(lsm_degenerate(__nwlsm_y, `__dir', `__fam')))
	if `"`seed'"' == "" local __seed = floor(runiform() * 2147483646) + 1
	else local __seed = `seed'
	tempname B V T D Z ACC PP MB ZM ME VT ES ER EO RH QM CM LM CLS WM PRM
	if !`__io' mata: __nwlsm_in = J(0, 0, .)
	capture noisily mata: LsmFit(__nwlsm_y, `__dir', __nwlsm_X, `dim', `nsample', `burnin', `interval', `__seed', `pilot', `plen', ///
		`=("`randomstart'" != "")', strtoreal(tokens("`__opt'")), ("`B'", "`V'", "`T'", "`D'", "`Z'", "`ACC'", "`PP'", "`MB'", "`ZM'", "`ME'", "`VT'", "`ES'", "`ER'", "`EO'", "`RH'", "`QM'", "`CM'", "`LM'", "`CLS'", "`WM'", "`PRM'"), __nwlsm_tr, __nwlsm_in)
	local __rc = _rc
	mata: st_local("__ny", strofreal(sum((__nwlsm_y :< .) :& (__nwlsm_y :!= .)) - sum(diagonal(__nwlsm_y) :< .)))
	capture mata: mata drop __nwlsm_y __nwlsm_X __nwlsm_tr __nwlsm_in
	if `__rc' exit `__rc'
	if !`__dir' local __ny = `__ny' / 2
	local __allnames `__names'
	local __eqs
	if `__io' {
		local __allnames `__names' `__onames'
		foreach __c of local __names {
			local __eqs `__eqs' network
		}
		foreach __c of local __onames {
			local __eqs `__eqs' outcome
		}
	}
	matrix colnames `B' = `__allnames'
	matrix rownames `V' = `__allnames'
	matrix colnames `V' = `__allnames'
	matrix rownames `T' = `__allnames'
	if `__io' {
		matrix coleq `B' = `__eqs'
		matrix roweq `V' = `__eqs'
		matrix coleq `V' = `__eqs'
		matrix roweq `T' = `__eqs'
	}
	matrix colnames `T' = mean sd mcse lower95 upper95
	local __effs
	if "`sender'" != "" local __effs `__effs' sender
	if "`receiver'" != "" local __effs `__effs' receiver
	if "`sociality'" != "" local __effs `__effs' sociality
	matrix colnames `ACC' = latent `__names' sender receiver sociality
	matrix rownames `ACC' = acceptance stepsize
	local __vn sigma2
	if `clusters' > 0 {
		local __vn
		forvalues __g = 1/`clusters' {
			local __vn `__vn' sigma2_`__g'
		}
	}
	foreach __e of local __effs {
		local __vn `__vn' `__e'_var
	}
	if `__fam' == 3 local __vn `__vn' dispersion
	if `__io' local __vn `__vn' sigma2_y
	matrix rownames `VT' = `__vn'
	matrix colnames `VT' = mean sd mcse lower95 upper95
	ereturn post `B' `V', obs(`__ny')
	ereturn local cmd "nwlsm"
	ereturn local version "1.0.5"
	ereturn local estat_cmd "nwlsm_estat"
	ereturn local title "`=cond(`__io', "Latent space model for social influence", "Latent space model")' (`=cond("`bilinear'" != "", "bilinear", "Euclidean distance")'; `__famname' family), posterior by MCMC"
	ereturn local latent = cond("`bilinear'" != "", "bilinear", "euclidean")
	ereturn local effects "`__effs'"
	ereturn local family "`__famname'"
	ereturn local network "`netname'"
	ereturn local directed = cond(`__dir', "directed", "undirected")
	ereturn local covariates "`__names'"
	if `__io' {
		ereturn local influence "yes"
		ereturn local outcome "`__ov'"
		ereturn local lagvar "`__lv'"
		ereturn local ocovariates "`__onames'"
		ereturn local weight "`__wtname'"
		ereturn scalar neighbors = cond(`neighbors' <= 0 | `neighbors' >= `nodes' - 1, `nodes' - 1, `neighbors')
		ereturn scalar obvar = `obvar'
		ereturn scalar osigma_a = `__osa'
		ereturn scalar osigma_b = `__osb'
		ereturn scalar prob_positive = `__pb2'
		ereturn scalar ll_outcome = `__ollk'
		ereturn scalar ll_outcome_sd = `__ollksd'
	}
	ereturn local positions "minimum Kullback-Leibler (MKL) positions"
	ereturn local start = cond("`randomstart'" != "", "random", cond(`__fam' == 0, "classical scaling of the geodesic distances", "classical scaling of the tie values"))
	ereturn scalar N_nodes = `nodes'
	ereturn scalar dim = `dim'
	ereturn scalar degenerate = `__degen'
	ereturn scalar nsample = `__nsample_total'
	ereturn scalar chains = `chains'
	ereturn scalar nsample_chain = `__nsample_total' / `chains'
	ereturn scalar threads = `=cond(`threads' > 0, `threads', 0)'
	ereturn scalar clusters = `clusters'
	if `__fam' == 2 & `"`trials'"' != "" {
		capture confirm number `trials'
		if !_rc ereturn scalar trials = `trials'
	}
	if `__fam' == 3 {
		ereturn scalar dispvar = `dispvar'
		ereturn scalar dispdf = `dispdf'
	}
	ereturn scalar burnin = `burnin'
	ereturn scalar interval = `interval'
	ereturn scalar seed = `__seed'
	ereturn scalar mklstarts = `mklstarts'
	ereturn scalar pilot = `pilot'
	ereturn scalar plen = `plen'
	ereturn scalar ll_mean = `__llk'
	ereturn scalar ll_sd = `__llksd'
	if `clusters' == 0 {
		ereturn scalar sigma2_mean = `__s2'
		ereturn scalar sigma2_sd = `__s2sd'
	}
	ereturn scalar prior_scale = `__priorscale'
	ereturn scalar prior_df = `__priordf'
	ereturn scalar mkl_start = `__mklstart'
	ereturn scalar mkl_iter = `__mkliter'
	ereturn scalar mkl_obj = `__mklobj'
	ereturn matrix summary = `T'
	ereturn matrix vsummary = `VT'
	ereturn matrix tuning = `ACC'
	if `chains' > 1 {
		local __rn `__names'
		foreach __c of local __onames {
			local __rn `__rn' y_`=cond("`__c'" == "_cons", "cons", "`__c'")'
		}
		local __rn `__rn' loglik `__vn'
		matrix rownames `RH' = `__rn'
		matrix colnames `RH' = rhat
		ereturn matrix rhat = `RH'
	}
	if `clusters' > 0 {
		local __dn
		forvalues __j = 1/`dim' {
			local __dn `__dn' dim`__j'
		}
		matrix colnames `CM' = `__dn'
		matrix colnames `LM' = mean sd
		ereturn matrix cmeans = `CM'
		ereturn matrix lambda = `LM'
		if `__matok' {
			local labs2 : subinstr local labs "," " ", all
			if `"`labs2'"' != "" capture matrix rownames `QM' = `labs2'
			if `"`labs2'"' != "" capture matrix rownames `CLS' = `labs2'
			ereturn matrix Q = `QM'
			matrix colnames `CLS' = class
			ereturn matrix class = `CLS'
		}
	}
	matrix colnames `MB' = `__names'
	matrix rownames `MB' = mkl
	ereturn matrix mkl_b = `MB'
	if `__matok' {
		local __dn
		forvalues __j = 1/`dim' {
			local __dn `__dn' dim`__j'
		}
		matrix colnames `Z' = `__dn'
		matrix colnames `ZM' = `__dn'
		local labs : subinstr local labs "," " ", all
		if `"`labs'"' != "" {
			foreach __m in Z ZM ME ES ER EO {
				capture matrix rownames ``__m'' = `labs'
			}
		}
		capture matrix rownames `PP' = `labs'
		capture matrix colnames `PP' = `labs'
		capture matrix rownames `D' = `labs'
		capture matrix colnames `D' = `labs'
		ereturn matrix D = `D'
		if `__io' {
			capture matrix rownames `WM' = `labs'
			capture matrix colnames `WM' = `labs'
			capture matrix rownames `PRM' = `labs'
			matrix colnames `PRM' = outcome predicted exposure influence
			ereturn matrix W = `WM'
			ereturn matrix pred = `PRM'
		}
		ereturn matrix Z = `Z'
		ereturn matrix Zmean = `ZM'
		ereturn matrix P = `PP'
		local __k 0
		foreach __e in sender receiver sociality {
			local ++__k
			if `: list __e in __effs' {
				if "`__e'" == "sender" {
					matrix colnames `ES' = mean sd
					ereturn matrix sender = `ES'
				}
				if "`__e'" == "receiver" {
					matrix colnames `ER' = mean sd
					ereturn matrix receiver = `ER'
				}
				if "`__e'" == "sociality" {
					matrix colnames `EO' = mean sd
					ereturn matrix sociality = `EO'
				}
			}
		}
		matrix colnames `ME' = sender receiver sociality
		ereturn matrix mkl_effects = `ME'
	}
	if `"`savedraws'"' != "" {
		preserve
		clear
		mata: LsmSaveDraws(tokens("`__names'"), strtoreal(tokens("`=("`sender'" != "")' `=("`receiver'" != "")' `=("`sociality'" != "")'")), tokens("`=subinstr("`__onames'", "_cons", "cons", 1)'"))
		save `"`savedraws'"', replace
		restore
	}
	_nwlsm_display
end

capture program drop _nwlsm_display
program _nwlsm_display
	di as text "{hline 78}"
	di as text "`e(title)'"
	di as text "Network: " as result "`e(network)'" as text " (" as result e(N_nodes) as text " actors, `e(directed)')" _col(52) as text "Latent dimensions: " as result e(dim)
	di as text "Draws: " as result e(nsample) as text "  burn-in: " as result e(burnin) as text "  interval: " as result e(interval) _col(52) as text "Seed: " as result e(seed)
	if e(chains) > 1 di as text "Chains: " as result e(chains) as text " (" as result e(nsample_chain) as text " draws each, pooled)" _c
	if e(clusters) > 0 {
		if e(chains) > 1 di as text _col(52) _c
		di as text "Latent clusters: " as result e(clusters)
	}
	else if e(chains) > 1 di
	local fmll = cond(abs(e(ll_mean)) < 10000, "%9.3f", "%10.2f")  // no scientific notation for large log likelihoods
	local fmlo = cond(abs(e(ll_outcome)) < 10000 | e(ll_outcome) >= ., "%9.3f", "%10.2f")
	di as text "Posterior mean log likelihood: " as result `fmll' e(ll_mean) as text " (sd " as result %6.3f e(ll_sd) as text ")" _c
	if "`e(influence)'" == "yes" di as text "  outcome: " as result `fmlo' e(ll_outcome) _c
	di
	di as text "{hline 78}"
	di as text _col(25) %10s "Mean" %10s "Sd" %10s "MC s.e." %11s "Lower 95%" %10s "Upper 95%"
	di as text "{hline 78}"
	tempname T
	matrix `T' = e(summary)
	local nm : rownames `T'
	local nnet = cond("`e(influence)'" == "yes", wordcount("`e(covariates)'"), .)
	local k 0
	if "`e(influence)'" == "yes" di as text "Network (ties)"
	foreach c of local nm {
		local ++k
		if `k' == `nnet' + 1 {
			di as text "{hline 78}"
			di as text "Outcome: " as result "`e(outcome)'" as text " (weights: `e(weight)', " as result e(neighbors) as text " nearest actors)"
		}
		di as text %-24s abbrev("`c'", 23) _col(25) as result %10.4f `T'[`k', 1] %10.4f `T'[`k', 2] %10.4f `T'[`k', 3] %11.4f `T'[`k', 4] %10.4f `T'[`k', 5]
	}
	di as text "{hline 78}"
	di as text "Variance parameters (posterior):"
	tempname VT
	matrix `VT' = e(vsummary)
	local vnm : rownames `VT'
	local k 0
	foreach c of local vnm {
		local ++k
		di as text %-24s "`c'" _col(25) as result %10.4f `VT'[`k', 1] %10.4f `VT'[`k', 2] %10.4f `VT'[`k', 3] %11.4f `VT'[`k', 4] %10.4f `VT'[`k', 5]
	}
	di as text "{hline 78}"
	tempname A
	matrix `A' = e(tuning)
	di as text "Acceptance rate of the latent positions: " as result %5.3f `A'[1, 1]
	if e(chains) > 1 {
		tempname RH
		matrix `RH' = e(rhat)
		mata: st_numscalar("__nwlsm_rmax", max(st_matrix("`RH'")))
		di as text "Gelman-Rubin statistic: largest value " as result %6.3f scalar(__nwlsm_rmax) as text " (see {bf:estat mcmcdiag})"
		capture scalar drop __nwlsm_rmax
	}
	if e(degenerate) == 1 di as error "{p 0 0 2}Warning: the network has no ties (or every possible tie); the data carry no information on the intercept, and its estimate is determined by the prior.{p_end}"
	if e(clusters) > 0 di as text "{p 0 0 2}Latent clusters: {bf:estat clusters} lists the posterior membership probabilities and the classification; cluster means {bf:e(cmeans)}, mixing proportions {bf:e(lambda)}.{p_end}"
	if "`e(influence)'" == "yes" di as text "{p 0 0 2}Influence model: {bf:estat influence} gives the influence effect, predictions and the actors who exert most influence; weights {bf:e(W)}.{p_end}"
	di as text "{p 0 0 2}Positions: minimum Kullback-Leibler estimate {bf:e(Z)} (also {bf:estat positions}); mean of the draws aligned to it {bf:e(Zmean)}; posterior mean distances {bf:e(D)}; posterior mean tie probabilities {bf:e(P)}. Positions are identified up to rotation, reflection and translation only.{p_end}"
end

capture mata: mata drop lsm_degenerate()
mata:
// 1 if the observed dyads (diagonal and, if undirected, the lower triangle left out) hold no tie or only ties (binary), or only zeros (poisson, binomial); else 0
real scalar lsm_degenerate(real matrix y, real scalar dir, real scalar fam)
{
	real scalar n, i, j, nobs, nz
	n = rows(y); nobs = 0; nz = 0
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j | (!dir & j < i) | y[i, j] >= .) continue
		nobs++
		nz = nz + (y[i, j] != 0)
	}
	if (nobs == 0 | fam == 3) return(0)
	if (nz == 0) return(1)
	if (fam == 0 & nz == nobs) return(1)
	return(0)
}
end
