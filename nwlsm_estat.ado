*! version 0.3.0  19nov2026 nwlsm_estat
* Copyright (C) 2026 Thomas U. Grund. Part of nwlsm; GNU GPL-3 (see LICENSE).

capture program drop nwlsm_estat
program nwlsm_estat, rclass
	version 14
	_nwlsm_libcheck
	gettoken subcmd rest : 0, parse(" ,")
	local subcmd = lower(trim("`subcmd'"))
	if "`e(cmd)'" != "nwlsm" error 301
	if "`subcmd'" == "positions" {
		_nwlsm_positions `rest'
		return add
	}
	else if "`subcmd'" == "mcmcdiag" {
		_nwlsm_mcmcdiag `rest'
		return add
	}
	else if "`subcmd'" == "clusters" {
		_nwlsm_clusters `rest'
		return add
	}
	else if "`subcmd'" == "predict" {
		_nwlsm_predict `rest'
		return add
	}
	else {
		di as err "estat `subcmd' not allowed after nwlsm"
		exit 321
	}
end

* estat positions: the positions (minimum Kullback-Leibler estimate; option mean: mean of the aligned draws), with the effects of the actors
capture program drop _nwlsm_positions
program _nwlsm_positions, rclass
	version 14
	syntax [, MEAN GENerate(string) NOList]
	tempname Z E
	if "`mean'" != "" matrix `Z' = e(Zmean)
	else matrix `Z' = e(Z)
	local n = rowsof(`Z')
	local d = colsof(`Z')
	local nm : rownames `Z'
	local cols
	local effs `"`e(effects)'"'
	local head "dim1"
	forvalues j = 2/`d' {
		local head `head' dim`j'
	}
	if `"`generate'"' != "" {
		if _N != `n' {
			di as error "{p 0 0 2}generate() needs a data set with one observation per actor in the order of the actors of {bf:`e(network)'}; the data set has `=_N' observations, the network `n' actors.{p_end}"
			exit 198
		}
		forvalues j = 1/`d' {
			confirm new variable `generate'`j'
		}
		forvalues j = 1/`d' {
			quietly generate double `generate'`j' = `Z'[_n, `j']
		}
		di as text "Positions stored in the variables {bf:`generate'1}" _c
		if `d' > 1 di as text " to {bf:`generate'`d'}" _c
		di as text "."
	}
	if "`nolist'" == "" {
		di
		di as text "Positions of the actors (" cond("`mean'" != "", "mean of the draws aligned to the minimum Kullback-Leibler estimate", "minimum Kullback-Leibler estimate") ")"
		local w = 12
		foreach c of local nm {
			local w = max(`w', min(20, udstrlen("`c'") + 1))
		}
		local line = `w' + 12 * `d' + 12 * wordcount("`effs'") + 2
		di as text "{hline `line'}"
		di as text %-`w's "Actor" _continue
		forvalues j = 1/`d' {
			di as text %12s "dim`j'" _continue
		}
		foreach e of local effs {
			di as text %12s "`e'" _continue
		}
		di
		di as text "{hline `line'}"
		forvalues i = 1/`n' {
			local lab : word `i' of `nm'
			if "`lab'" == "" local lab "`i'"
			di as text %-`w's abbrev("`lab'", `w' - 1) _continue
			forvalues j = 1/`d' {
				di as result %12.4f `Z'[`i', `j'] _continue
			}
			foreach e of local effs {
				matrix `E' = e(`e')
				di as result %12.4f `E'[`i', 1] _continue
			}
			di
		}
		di as text "{hline `line'}"
		if "`effs'" != "" di as text "{p 0 0 2}Effects: posterior means. Positions are identified up to rotation, reflection and translation only.{p_end}"
		else di as text "{p 0 0 2}Positions are identified up to rotation, reflection and translation only.{p_end}"
	}
	return matrix Z = `Z'
end

* estat mcmcdiag: effective sample size and Geweke z-score of the coefficients, the log likelihood and the variance parameters
capture program drop _nwlsm_mcmcdiag
program _nwlsm_mcmcdiag, rclass
	version 14
	syntax [, NOList]
	capture mata: st_numscalar("__nwlsm_Sck", rows(LsmSample("beta")))
	if _rc {
		di as error "the posterior draws of the last fit are no longer in memory; fit the model again"
		exit 498
	}
	scalar drop __nwlsm_Sck
	tempname R
	local C = e(chains)
	if `C' >= . local C 1
	local G = e(clusters)
	if `G' >= . local G 0
	mata: __nwlsm_dg = LsmDiag((LsmSample("beta"), LsmSample("llk"), LsmSample("s2")), `C')
	local nm : colnames e(b)
	if `G' > 0 {
		local nm `nm' loglik
		forvalues g = 1/`G' {
			local nm `nm' sigma2_`g'
		}
	}
	else local nm `nm' loglik sigma2
	local act = 0
	local effs `"`e(effects)'"'
	foreach e of local effs {
		local ++act
	}
	if `act' > 0 {
		local cols
		local k 0
		foreach t in sender receiver sociality {
			local ++k
			if `: list t in effs' {
				mata: __nwlsm_dg = __nwlsm_dg \ LsmDiag(LsmSample("ev")[., `k'], `C')
				local nm `nm' `t'_var
			}
		}
	}
	mata: st_matrix("`R'", __nwlsm_dg)
	capture mata: mata drop __nwlsm_dg
	matrix rownames `R' = `nm'
	matrix colnames `R' = ess geweke mean sd rhat
	local S = e(nsample)
	if "`nolist'" == "" {
		di
		di as text "MCMC diagnostics (" as result e(nsample) as text " draws" _c
		if `C' > 1 di as text " in " as result `C' as text " chains" _c
		di as text ", interval " as result e(interval) as text ")"
		di as text "{hline 76}"
		di as text %-22s "Parameter" %10s "Mean" %10s "Sd" %10s "ESS" %10s "ESS/draw" %10s "Geweke z" _c
		if `C' > 1 di as text %10s "Rhat" _c
		di
		di as text "{hline 76}"
		local k 0
		foreach c of local nm {
			local ++k
			di as text %-22s abbrev("`c'", 21) as result %10.4f `R'[`k', 3] %10.4f `R'[`k', 4] %10.1f `R'[`k', 1] %10.3f `R'[`k', 1] / `S' %10.2f `R'[`k', 2] _c
			if `C' > 1 di as result %10.3f `R'[`k', 5] _c
			di
		}
		di as text "{hline 76}"
		tempname A
		matrix `A' = e(tuning)
		di as text "Acceptance rates: positions " as result %5.3f `A'[1, 1] _c
		local k 1
		foreach c in `: colnames e(b)' {
			local ++k
			di as text ", `c' " as result %5.3f `A'[1, `k'] _c
		}
		local k = colsof(e(b)) + 1
		local j 0
		foreach t in sender receiver sociality {
			local ++j
			if `: list t in effs' di as text ", `t' " as result %5.3f `A'[1, `k' + `j'] _c
		}
		di
		di as text "{p 0 0 2}ESS: effective sample size by the initial monotone positive sequence estimator. Geweke z: first 10% of the draws against the last 50%; |z| above 2 suggests that the chain has not settled.{p_end}"
		if `C' > 1 di as text "{p 0 0 2}With `C' chains: ESS is the sum over the chains, Geweke z the value of largest size over the chains, Rhat the Gelman-Rubin potential scale reduction factor (values close to 1, say below 1.1, indicate that the chains agree).{p_end}"
	}
	return matrix diag = `R'
end

* estat predict: the posterior mean tie probabilities as a valued network of nwcommands (and r(P))
capture program drop _nwlsm_predict
program _nwlsm_predict, rclass
	version 14
	syntax [, GENerate(string) REPLACE]
	tempname P
	matrix `P' = e(P)
	local n = rowsof(`P')
	if `"`generate'"' != "" {
		capture _nwsyntax `generate', max(1)
		local exists = (_rc == 0)
		if `exists' & "`replace'" == "" {
			di as error "{p 0 0 2}network {bf:`generate'} exists already; use another name or the replace option.{p_end}"
			exit 110
		}
		if `exists' nwdrop `generate'
		preserve
		clear
		quietly set obs `n'
		forvalues j = 1/`n' {
			quietly gen double __lsmp`j' = .
		}
		mata: st_store(., ., st_matrix("`P'"))
		mata: st_matrix("__lsmP", st_matrix("`P'"))
		if "`e(directed)'" == "directed" nwset __lsmp*, name(`generate') directed valued
		else nwset __lsmp*, name(`generate') undirected valued
		restore
		capture matrix drop __lsmP
		di as text "Posterior mean tie probabilities stored as the valued network {bf:`generate'}."
	}
	return matrix P = `P'
end

* estat clusters: posterior membership probabilities and the hard classification (most probable cluster) of the latent cluster model
capture program drop _nwlsm_clusters
program _nwlsm_clusters, rclass
	version 14
	syntax [, GENerate(string) NOList]
	if "`e(clusters)'" == "" | "`e(clusters)'" == "0" {
		di as error "the model was fitted without latent clusters (option clusters())"
		exit 198
	}
	tempname Q CL
	capture matrix `Q' = e(Q)
	if _rc {
		di as error "the membership matrix was not stored (too many actors for a Stata matrix); see LsmSample(""class"")"
		exit 498
	}
	matrix `CL' = e(class)
	local n = rowsof(`Q')
	local G = colsof(`Q')
	local nm : rownames `Q'
	if `"`generate'"' != "" {
		if _N != `n' {
			di as error "{p 0 0 2}generate() needs a data set with one observation per actor in the order of the actors of {bf:`e(network)'}; the data set has `=_N' observations, the network `n' actors.{p_end}"
			exit 198
		}
		confirm new variable `generate'_class
		forvalues g = 1/`G' {
			confirm new variable `generate'_p`g'
		}
		quietly generate byte `generate'_class = `CL'[_n, 1]
		forvalues g = 1/`G' {
			quietly generate double `generate'_p`g' = `Q'[_n, `g']
		}
		di as text "Classification stored in {bf:`generate'_class}, probabilities in {bf:`generate'_p1} to {bf:`generate'_p`G'}."
	}
	if "`nolist'" == "" {
		di
		di as text "Latent clusters: posterior membership probabilities (labels corrected over the draws)"
		local line = 16 + 8 + 10 * `G'
		di as text "{hline `line'}"
		di as text %-16s "Actor" %8s "Class" _c
		forvalues g = 1/`G' {
			di as text %10s "P(`g')" _c
		}
		di
		di as text "{hline `line'}"
		forvalues i = 1/`n' {
			local lab : word `i' of `nm'
			if "`lab'" == "" local lab "`i'"
			di as text %-16s abbrev("`lab'", 15) as result %8.0f `CL'[`i', 1] _c
			forvalues g = 1/`G' {
				di as result %10.3f `Q'[`i', `g'] _c
			}
			di
		}
		di as text "{hline `line'}"
		di as text "{p 0 0 2}Class: the cluster with the largest posterior probability.{p_end}"
	}
	return matrix Q = `Q'
	return matrix class = `CL'
end
