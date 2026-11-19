*! version 0.3.0  19nov2026 nwlsm
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
		SENDER RECEIVER SOCiality SENDERVar(real 1) SENDERDf(real 3) RECEIVERVar(real 1) RECEIVERDf(real 3) SOCIALITYVar(real 1) SOCIALITYDf(real 3) BILinear MKLStarts(integer 1) CHAins(integer 1) THReads(integer 0) CLUsters(integer 0)]
	if `dim' < 1 | `dim' > 10 {
		di as error "dim() must be an integer between 1 and 10"
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
	if `mklstarts' < 1 {
		di as error "mklstarts() must be at least 1"
		error 198
	}
	if `chains' < 1 | `chains' > 64 {
		di as error "chains() must be an integer between 1 and 64"
		error 198
	}
	if `threads' < 0 | `threads' > 256 {
		di as error "threads() must be between 1 and 256 (default: one thread per core, at most one per chain)"
		error 198
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
	_nwlsm_libcheck
	set more off
	_nwsyntax `netname', max(1)
	if "`is2mode'" == "true" {
		di as error "nwlsm fits one-mode networks; `netname' is a two-mode network."
		error 198
	}
	if "`valued'" == "true" {
		di as error "{p 0 0 2}nwlsm fits binary networks; `netname' is a valued network. Use {bf:nwset} or {bf:nwreplace} to dichotomise it first.{p_end}"
		error 198
	}
	if `nodes' < 3 {
		di as error "the network needs at least 3 actors"
		error 198
	}
	local __dir = ("`directed'" == "true")
	if !`__dir' & ("`sender'" != "" | "`receiver'" != "") {
		di as error "{p 0 0 2}sender and receiver effects need a directed network; use {bf:sociality} for the undirected network {bf:`netname'}.{p_end}"
		error 198
	}
	if "`sociality'" != "" & ("`sender'" != "" | "`receiver'" != "") & !`__dir' {
		di as error "sociality cannot be combined with sender or receiver effects in an undirected network"
		error 198
	}
	local __opt "`=("`sender'" != "")' `=("`receiver'" != "")' `=("`sociality'" != "")' `=("`bilinear'" != "")' `sendervar' `senderdf' `receivervar' `receiverdf' `socialityvar' `socialitydf' `mklstarts' `chains' `threads' `clusters'"
	capture mata: mata drop __nwlsm_y __nwlsm_X
	mata: __nwlsm_y = lsm_ymat(*(`netobj'->get_matrix()))
	mata: __nwlsm_X = J(0, `nodes', .)
	local __names "_cons"
	foreach __c of local edgecov {
		capture _nwsyntax `__c', max(1) other(__cv)
		if !_rc {
			if "`__cvis2mode'" == "true" | `__cvnodes' != `nodes' {
				capture mata: mata drop __nwlsm_y __nwlsm_X
				di as error "{p 0 0 2}the edge covariate {bf:`__c'} must be a one-mode network of the same `nodes' actors as {bf:`netname'}.{p_end}"
				error 198
			}
			mata: __nwlsm_X = __nwlsm_X \ *(`__cvnetobj'->get_matrix())
		}
		else {
			capture confirm matrix `__c'
			if _rc {
				capture mata: mata drop __nwlsm_y __nwlsm_X
				di as error "{p 0 0 2}edgecov(): {bf:`__c'} is neither a network nor a Stata matrix.{p_end}"
				error 198
			}
			if rowsof(`__c') != `nodes' | colsof(`__c') != `nodes' {
				capture mata: mata drop __nwlsm_y __nwlsm_X
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
			capture mata: mata drop __nwlsm_y __nwlsm_X
			di as error "absdiff(): not a list of variables in the current dataset"
			error 111
		}
		capture confirm numeric variable `absdiff'
		if _rc {
			capture mata: mata drop __nwlsm_y __nwlsm_X
			di as error "absdiff(): the variables must be numeric"
			error 109
		}
		if _N != `nodes' {
			capture mata: mata drop __nwlsm_y __nwlsm_X
			di as error "{p 0 0 2}absdiff() needs one observation per actor in the current dataset, in the order of the actors of {bf:`netname'} (the dataset has `=_N' observations, the network `nodes' actors).{p_end}"
			error 198
		}
		foreach __v of local absdiff {
			mata: st_local("__hm", strofreal(hasmissing(st_data(., "`__v'"))))
			if `__hm' {
				capture mata: mata drop __nwlsm_y __nwlsm_X
				di as error "absdiff(): variable {bf:`__v'} has missing values"
				error 416
			}
			mata: __nwlsm_X = __nwlsm_X \ lsm_node_covs("`__v'", `nodes')
			local __names `__names' absdiff_`__v'
		}
	}
	mata: st_local("__hm", strofreal(hasmissing(__nwlsm_X)))
	if `__hm' {
		capture mata: mata drop __nwlsm_y __nwlsm_X
		di as error "edge covariates must not have missing values"
		error 416
	}
	local __nc : word count `__names'
	local __dp : list dups __names
	if "`__dp'" != "" {
		capture mata: mata drop __nwlsm_y __nwlsm_X
		di as error "the same edge covariate is listed twice"
		error 198
	}
	if `"`seed'"' == "" local __seed = floor(runiform() * 2147483646) + 1
	else local __seed = `seed'
	tempname B V T D Z ACC PP MB ZM ME VT ES ER EO RH QM CM LM CLS
	capture noisily mata: LsmFit(__nwlsm_y, `__dir', __nwlsm_X, `dim', `nsample', `burnin', `interval', `__seed', `pilot', `plen', ///
		`=("`randomstart'" != "")', strtoreal(tokens("`__opt'")), ("`B'", "`V'", "`T'", "`D'", "`Z'", "`ACC'", "`PP'", "`MB'", "`ZM'", "`ME'", "`VT'", "`ES'", "`ER'", "`EO'", "`RH'", "`QM'", "`CM'", "`LM'", "`CLS'"))
	local __rc = _rc
	mata: st_local("__ny", strofreal(sum((__nwlsm_y :< .) :& (__nwlsm_y :!= .)) - sum(diagonal(__nwlsm_y) :< .)))
	capture mata: mata drop __nwlsm_y __nwlsm_X
	if `__rc' exit `__rc'
	if !`__dir' local __ny = `__ny' / 2
	matrix colnames `B' = `__names'
	matrix rownames `V' = `__names'
	matrix colnames `V' = `__names'
	matrix rownames `T' = `__names'
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
	matrix rownames `VT' = `__vn'
	matrix colnames `VT' = mean sd mcse lower95 upper95
	ereturn post `B' `V', obs(`__ny')
	ereturn local cmd "nwlsm"
	ereturn local version "0.4.0"
	ereturn local estat_cmd "nwlsm_estat"
	ereturn local title "Latent space model (`=cond("`bilinear'" != "", "bilinear", "Euclidean distance")'), posterior by MCMC"
	ereturn local latent = cond("`bilinear'" != "", "bilinear", "euclidean")
	ereturn local effects "`__effs'"
	ereturn local network "`netname'"
	ereturn local directed = cond(`__dir', "directed", "undirected")
	ereturn local covariates "`__names'"
	ereturn local positions "minimum Kullback-Leibler (MKL) positions"
	ereturn local start = cond("`randomstart'" != "", "random", "classical scaling of the geodesic distances")
	ereturn scalar N_nodes = `nodes'
	ereturn scalar dim = `dim'
	ereturn scalar nsample = `__nsample_total'
	ereturn scalar chains = `chains'
	ereturn scalar nsample_chain = `__nsample_total' / `chains'
	ereturn scalar threads = `=cond(`threads' > 0, `threads', 0)'
	ereturn scalar clusters = `clusters'
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
		local __rn : rownames e(summary)
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
		mata: LsmSaveDraws(tokens("`__names'"), strtoreal(tokens("`=("`sender'" != "")' `=("`receiver'" != "")' `=("`sociality'" != "")'")))
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
	di as text "Posterior mean log likelihood: " as result %9.3f e(ll_mean) as text " (sd " as result %6.3f e(ll_sd) as text ")"
	di as text "{hline 78}"
	di as text _col(25) %10s "Mean" %10s "Sd" %10s "MC s.e." %11s "Lower 95%" %10s "Upper 95%"
	di as text "{hline 78}"
	tempname T
	matrix `T' = e(summary)
	local nm : rownames `T'
	local k 0
	foreach c of local nm {
		local ++k
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
	if e(clusters) > 0 di as text "{p 0 0 2}Latent clusters: {bf:estat clusters} lists the posterior membership probabilities and the classification; cluster means {bf:e(cmeans)}, mixing proportions {bf:e(lambda)}.{p_end}"
	di as text "{p 0 0 2}Positions: minimum Kullback-Leibler estimate {bf:e(Z)} (also {bf:estat positions}); mean of the draws aligned to it {bf:e(Zmean)}; posterior mean distances {bf:e(D)}; posterior mean tie probabilities {bf:e(P)}. Positions are identified up to rotation, reflection and translation only.{p_end}"
end
