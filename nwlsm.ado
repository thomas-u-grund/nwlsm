*! version 0.2.0  19nov2026 nwlsm
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
		INTerval(integer 10) SEED(string) PILOT(integer 4) PLEN(integer 500) RANDOMstart SAVEDraws(string)]
	if `dim' < 1 | `dim' > 10 {
		di as error "dim() must be an integer between 1 and 10"
		error 198
	}
	if `nsample' < 10 | `burnin' < 0 | `interval' < 1 | `pilot' < 0 | `plen' < 10 {
		di as error "{p 0 0 2}nsample() must be at least 10, burnin() not negative, interval() at least 1, pilot() not negative and plen() at least 10.{p_end}"
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
	tempname B V T D Z ACC
	capture noisily mata: LsmFit(__nwlsm_y, `__dir', __nwlsm_X, `dim', `nsample', `burnin', `interval', `__seed', `pilot', `plen', ///
		`=("`randomstart'" != "")', "`B'", "`V'", "`T'", "`D'", "`Z'", "`ACC'")
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
	local __ac `__names' 
	matrix colnames `ACC' = latent `__names'
	matrix rownames `ACC' = acceptance stepsize
	ereturn post `B' `V', obs(`__ny')
	ereturn local cmd "nwlsm"
	ereturn local version "0.2.0"
	ereturn local title "Latent space model (Euclidean distance), posterior by MCMC"
	ereturn local network "`netname'"
	ereturn local directed = cond(`__dir', "directed", "undirected")
	ereturn local covariates "`__names'"
	ereturn local positions "mean of the Procrustes-aligned draws"
	ereturn local start = cond("`randomstart'" != "", "random", "classical scaling of the geodesic distances")
	ereturn scalar N_nodes = `nodes'
	ereturn scalar dim = `dim'
	ereturn scalar nsample = `nsample'
	ereturn scalar burnin = `burnin'
	ereturn scalar interval = `interval'
	ereturn scalar seed = `__seed'
	ereturn scalar pilot = `pilot'
	ereturn scalar plen = `plen'
	ereturn scalar ll_mean = `__llk'
	ereturn scalar ll_sd = `__llksd'
	ereturn scalar sigma2_mean = `__s2'
	ereturn scalar sigma2_sd = `__s2sd'
	ereturn scalar prior_scale = `__priorscale'
	ereturn scalar prior_df = `__priordf'
	ereturn matrix summary = `T'
	ereturn matrix tuning = `ACC'
	if `__matok' {
		matrix colnames `Z' = dim1-dim`dim'
		ereturn matrix D = `D'
		ereturn matrix Z = `Z'
	}
	if `"`savedraws'"' != "" {
		preserve
		clear
		mata: LsmSaveDraws(tokens("`__names'"))
		qui rename _cons beta_cons
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
	di as text "Variance of the latent positions (sigma2), posterior mean: " as result %7.3f e(sigma2_mean)
	tempname A
	matrix `A' = e(tuning)
	di as text "Acceptance rate of the latent positions: " as result %5.3f `A'[1, 1]
	di as text "{p 0 0 2}Posterior mean distances: {bf:e(D)}; mean positions (aligned by Procrustes rotation): {bf:e(Z)}. Positions are identified up to rotation, reflection and translation.{p_end}"
end
