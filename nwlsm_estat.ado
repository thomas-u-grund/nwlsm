*! version 0.6.0  19nov2026 nwlsm_estat
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
	else if "`subcmd'" == "gof" {
		_nwlsm_gof `rest'
		return add
	}
	else if "`subcmd'" == "influence" {
		_nwlsm_influence `rest'
		return add
	}
	else if "`subcmd'" == "ppcheck" {
		_nwlsm_ppcheck `rest'
		return add
	}
	else if "`subcmd'" == "bic" {
		_nwlsm_bic `rest'
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
	syntax [, NOList GRAPH NAME(string)]
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
	mata: __nwlsm_dg = LsmDiag((LsmSample("beta"), LsmSample("ob"), LsmSample("llk"), LsmSample("s2")), `C')
	local nm `e(covariates)'
	if "`e(influence)'" == "yes" {
		foreach c in `e(ocovariates)' {
			local nm `nm' y_`=cond("`c'" == "_cons", "cons", "`c'")'
		}
	}
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
	if "`e(family)'" == "normal" {
		mata: __nwlsm_dg = __nwlsm_dg \ LsmDiag(LsmSample("disp"), `C')
		local nm `nm' dispersion
	}
	if "`e(influence)'" == "yes" {
		mata: __nwlsm_dg = __nwlsm_dg \ LsmDiag(LsmSample("os2"), `C')
		local nm `nm' sigma2_y
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
		foreach c in `e(covariates)' {
			local ++k
			di as text ", `c' " as result %5.3f `A'[1, `k'] _c
		}
		local k = wordcount("`e(covariates)'") + 1
		local j 0
		foreach t in sender receiver sociality {
			local ++j
			if `: list t in effs' di as text ", `t' " as result %5.3f `A'[1, `k' + `j'] _c
		}
		di
		di as text "{p 0 0 2}ESS: effective sample size by the initial monotone positive sequence estimator. Geweke z: first 10% of the draws against the last 50%; |z| above 2 suggests that the chain has not settled.{p_end}"
		if `C' > 1 di as text "{p 0 0 2}With `C' chains: ESS is the sum over the chains, Geweke z the value of largest size over the chains, Rhat the Gelman-Rubin potential scale reduction factor (values close to 1, say below 1.1, indicate that the chains agree).{p_end}"
	}
	if "`graph'" != "" {
		if "`name'" == "" local name "mcmcdiag"
		local np : word count `nm'
		local a1 = (`: list posof "sender" in effs' > 0)
		local a2 = (`: list posof "receiver" in effs' > 0)
		local a3 = (`: list posof "sociality" in effs' > 0)
		mata: st_matrix("__nwlsm_trace", LsmTraceMat((`a1', `a2', `a3')))
		local ttl "Trace plots of the posterior draws"
		if `C' > 1 local ttl "`ttl' (`C' chains)"
		local Sc = rowsof(__nwlsm_trace) / `C'
		preserve
		quietly drop _all
		quietly set obs `Sc'
		quietly gen draw = _n
		local panels
		forvalues k = 1/`np' {
			local pn : word `k' of `nm'
			local lines
			forvalues c = 1/`C' {
				quietly gen double v`k'_`c' = .
				mata: st_store(., "v`k'_`c'", st_matrix("__nwlsm_trace")[((`c' - 1) * `Sc' + 1)..(`c' * `Sc'), `k'])
				local col = word("navy maroon forest_green dkorange purple teal sienna gs8", mod(`c' - 1, 8) + 1)
				local lines `lines' (line v`k'_`c' draw, lcolor(`col') lwidth(vthin))
			}
			local sub
			if `C' > 1 local sub subtitle("Rhat = `: di %5.3f `R'[`k', 5]'", size(vsmall))
			twoway `lines', title("`pn'", size(small)) `sub' xtitle("draw") ytitle("") legend(off) nodraw name(__lsmtr`k', replace)
			local panels `panels' __lsmtr`k'
		}
		local cols = cond(`np' <= 4, 2, 3)
		graph combine `panels', cols(`cols') title("`ttl'", size(medium)) name(`name', replace)
		forvalues k = 1/`np' {
			capture graph drop __lsmtr`k'
		}
		restore
		capture matrix drop __nwlsm_trace
		di as txt "{p 0 0 2}(trace plots saved as {bf:`name'}; one line per chain)"
	}
	return matrix diag = `R'
end

* estat predict: the posterior mean tie probabilities (expected values of a valued family) as a valued network of nwcommands (and r(P))
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
		di as text cond("`e(family)'" == "poisson" | "`e(family)'" == "binomial" | "`e(family)'" == "normal", "Posterior mean expected values stored as the valued network {bf:`generate'}.", "Posterior mean tie probabilities stored as the valued network {bf:`generate'}.")
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

* estat gof: goodness of fit as latentnet's gof(): networks simulated from posterior draws against the observed network
capture program drop _nwlsm_gof
program _nwlsm_gof, rclass
	version 14
	syntax [, NSIM(integer 100) SEED(string) GOF(string) PLOT MAXDEG(integer -1) MAXDIST(integer -1) NAME(string) THReads(integer 0)]
	capture mata: st_numscalar("__nwlsm_Sck", rows(LsmSample("beta")))
	if _rc {
		di as error "the posterior draws of the last fit are no longer in memory; fit the model again"
		exit 498
	}
	scalar drop __nwlsm_Sck
	if `nsim' < 2 {
		di as error "nsim() must be at least 2"
		exit 198
	}
	if `threads' < 0 | `threads' > 256 {
		di as error "threads() must be between 1 and 256"
		exit 198
	}
	if `"`seed'"' != "" {
		capture confirm integer number `seed'
		if _rc | real(`"`seed'"') < 0 | real(`"`seed'"') > 2147483647 {
			di as error "seed() must be an integer between 0 and 2147483647"
			exit 198
		}
	}
	else local seed = floor(runiform() * 2147483646) + 1
	local dir = ("`e(directed)'" == "directed")
	local fam "`e(family)'"
	if "`fam'" == "" local fam "bernoulli"
	local valued = ("`fam'" != "bernoulli")
	local normal = ("`fam'" == "normal")
	* statistics: R's  GOF = ~ degree + ... ; default as latentnet (idegree + odegree + distance; degree + distance if undirected), plus the values of a valued network
	if `"`gof'"' == "" {
		if `normal' local gof "model values"
		else if `dir' local gof "idegree odegree distance"
		else local gof "degree distance"
		if `valued' & !`normal' local gof "`gof' values"
	}
	local gof = subinstr(subinstr(`"`gof'"', "~", " ", .), "+", " ", .)
	local allterms "degree idegree odegree espartners dspartners distance model values"
	local mask 0
	local want
	foreach tok of local gof {
		local found
		foreach t of local allterms {
			if "`tok'" == "`t'" local found "`t'"
		}
		if "`found'" == "" {
			local nm 0
			foreach t of local allterms {
				if strpos("`t'", "`tok'") == 1 {
					local found "`t'"
					local ++nm
				}
			}
			if `nm' != 1 local found
		}
		if "`found'" == "" {
			di as error "{p 0 0 2}estat gof: unknown statistic {bf:`tok'}; choose from degree, idegree, odegree, espartners, dspartners, distance, model and values.{p_end}"
			exit 198
		}
		if !`dir' & inlist("`found'", "idegree", "odegree") {
			di as error "estat gof: `found' needs a directed network; use degree."
			exit 198
		}
		if `dir' & "`found'" == "degree" {
			di as error "estat gof: degree is for undirected networks; use idegree and odegree."
			exit 198
		}
		if "`found'" == "values" & !`valued' {
			di as error "estat gof: values is for the valued families (family(poisson), family(binomial) or family(normal))."
			exit 198
		}
		if `normal' & !inlist("`found'", "model", "values") {
			di as error "{p 0 0 2}estat gof: after family(normal) the statistics of the network structure are not defined; use model and values.{p_end}"
			exit 198
		}
		local inw : list found in want
		if !`inw' local want `want' `found'
	}
	foreach t of local want {
		if "`t'" == "degree" | "`t'" == "idegree" local mask = `mask' + 1
		if "`t'" == "odegree" local mask = `mask' + 2
		if "`t'" == "espartners" local mask = `mask' + 4
		if "`t'" == "dspartners" local mask = `mask' + 8
		if "`t'" == "distance" local mask = `mask' + 16
		if "`t'" == "model" local mask = `mask' + 32
		if "`t'" == "values" local mask = `mask' + 64
	}
	mata: st_numscalar("__nwlsm_S", rows(LsmSample("beta")))
	if `nsim' > scalar(__nwlsm_S) {
		di as text "(nsim() reduced to the number of posterior draws, " scalar(__nwlsm_S) ")"
		local nsim = scalar(__nwlsm_S)
	}
	scalar drop __nwlsm_S
	mata: LsmGof(`nsim', `seed', `mask', `threads', 10, "__lsmgof")
	* print
	local order "degree idegree odegree espartners dspartners distance model values"
	local g_degree deg
	local g_idegree ideg
	local g_odegree odeg
	local g_espartners esp
	local g_dspartners dsp
	local g_distance dist
	local g_model model
	local g_values val
	local d_degree "degree"
	local d_idegree "in-degree"
	local d_odegree "out-degree"
	local d_espartners "edgewise shared partners"
	local d_dspartners "dyadwise shared partners"
	local d_distance "minimum geodesic distance"
	local d_model "model statistics"
	local d_values "tie values"
	local cnames `e(covariates)'
	foreach t of local order {
		local inw : list t in want
		if !`inw' continue
		local g "`g_`t''"
		tempname P
		matrix `P' = __lsmgof_`g'
		local k = rowsof(`P')
		local labs
		forvalues i = 1/`k' {
			if "`t'" == "model" {
				if `i' == 1 local lb = cond(`valued', "total_value", "edges")
				else {
					local lb : word `i' of `cnames'
				}
			}
			else if "`t'" == "distance" local lb = cond(`i' == `k', "Inf", "`i'")
			else if "`t'" == "values" {
				if `normal' local lb = "bin`i'"
				else local lb = cond(`i' == `k', "`=`k'-1'+", "`=`i'-1'")
			}
			else local lb = `i' - 1
			local labs `labs' `lb'
		}
		local hd "Statistic"
		if "`t'" == "degree" local hd "Degree"
		if "`t'" == "idegree" local hd "In-degree"
		if "`t'" == "odegree" local hd "Out-degree"
		if "`t'" == "espartners" local hd "ESP"
		if "`t'" == "dspartners" local hd "DSP"
		if "`t'" == "distance" local hd "Distance"
		if "`t'" == "values" local hd = cond(`normal', "Value bin", "Value")
		local w = 12
		foreach lb of local labs {
			local w = max(`w', min(28, udstrlen("`lb'")))
		}
		matrix colnames `P' = obs min mean max MC_p-value
		matrix rownames `P' = `labs'
		di
		di as txt "Goodness of fit: `d_`t'' (" as res `nsim' as txt " networks simulated from the posterior)"
		di as txt "{hline `=`w'+1'}{c TT}{hline 52}"
		di as txt %`w's "`hd'" " {c |}" %10s "Observed" %10s "Min" %10s "Mean" %10s "Max" %12s "MC p-value"
		di as txt "{hline `=`w'+1'}{c +}{hline 52}"
		forvalues i = 1/`k' {
			if `P'[`i',1] == 0 & `P'[`i',2] == 0 & `P'[`i',4] == 0 continue
			local lb : word `i' of `labs'
			local fm "%10.0f"
			if `P'[`i',1] != floor(`P'[`i',1]) | `P'[`i',2] != floor(`P'[`i',2]) | `P'[`i',4] != floor(`P'[`i',4]) local fm "%10.2f"
			di as txt %`w's abbrev("`lb'", `w') " {c |}" as res `fm' `P'[`i',1] `fm' `P'[`i',2] %10.2f `P'[`i',3] `fm' `P'[`i',4] %12.2f `P'[`i',5]
		}
		di as txt "{hline `=`w'+1'}{c BT}{hline 52}"
		return matrix gof_`t' = `P'
	}
	di as txt "{p 0 0 2}Each simulated network is drawn from the model at one of `nsim' posterior draws (evenly spaced over the stored draws). MC p-value: 2 min(P(simulated <= observed), P(simulated >= observed)), at most 1, as in latentnet. Missing dyads of the observed network count as no tie.{p_end}"
	if `valued' & !`normal' di as txt "{p 0 0 2}Degree, shared partners and distances are those of the network of the dyads with a positive value.{p_end}"
	if "`plot'" != "" {
		if "`name'" == "" local name "gof"
		local panels
		local np 0
		foreach t in degree idegree odegree distance values {
			local inw : list t in want
			if !`inw' continue
			local g "`g_`t''"
			local ++np
			local startval = cond("`t'" == "distance", 1, 0)
			local ncell = rowsof(__lsmgof_q_`g')
			* shown: the cells up to the last one with a positive observed or simulated share (at most maxdeg() + 1 or maxdist() cells if given)
			mata: st_local("maxc", strofreal(max(selectindex(rowmax((abs(st_matrix("__lsmgof_q_`g'")[., 1]), st_matrix("__lsmgof_q_`g'")[., 6])) :> 0) \ 1)))
			if "`t'" == "distance" & `maxdist' >= 0 local maxc = min(`maxc', `maxdist')
			if inlist("`t'", "degree", "idegree", "odegree") & `maxdeg' >= 0 local maxc = min(`maxc', `maxdeg' + 1)
			local lastlab = cond("`t'" == "distance", "NR", "`=`maxc'-1'")
			_nwlsm_gofplot __lsmgof_q_`g' "`d_`t''" "panel`np'" `startval' `maxc' "`t'" `ncell' "`lastlab'"
			local panels `panels' panel`np'
		}
		if `np' > 0 {
			graph combine `panels', cols(`=min(`np', 2)') title("Goodness of fit: observed against `nsim' simulated networks", size(medium)) name(`name', replace)
			foreach pn of local panels {
				capture graph drop `pn'
			}
			di as txt "{p 0 0 2}(plot saved as {bf:`name'}; boxes and whiskers summarize the simulated networks, the red line with circles is the observed network){p_end}"
		}
	}
	return scalar nsim = `nsim'
	return scalar seed = `seed'
	foreach m in deg ideg odeg esp dsp dist model val {
		capture matrix drop __lsmgof_`m'
		capture matrix drop __lsmgof_q_`m'
	}
	capture scalar drop __lsmgof_nsim
end

* one panel of the estat gof plot: observed share (red) and the box and whiskers of the simulated shares; matrix columns obs, min, p25, median, p75, max
capture program drop _nwlsm_gofplot
program _nwlsm_gofplot
	args mat ttl graphname startval nshow stat ncell lastlab
	preserve
	quietly drop _all
	quietly set obs `nshow'
	tempname Q
	matrix `Q' = `mat'
	mata: st_store(., st_addvar("double", "observed"), st_matrix("`Q'")[1..`nshow', 1])
	mata: st_store(., st_addvar("double", "ymin"), st_matrix("`Q'")[1..`nshow', 2])
	mata: st_store(., st_addvar("double", "yp25"), st_matrix("`Q'")[1..`nshow', 3])
	mata: st_store(., st_addvar("double", "ymedian"), st_matrix("`Q'")[1..`nshow', 4])
	mata: st_store(., st_addvar("double", "yp75"), st_matrix("`Q'")[1..`nshow', 5])
	mata: st_store(., st_addvar("double", "ymax"), st_matrix("`Q'")[1..`nshow', 6])
	quietly gen xval = `startval' + _n - 1
	local xlab
	forvalues i = 1/`nshow' {
		local v = `startval' + `i' - 1
		local lb "`v'"
		if "`stat'" == "distance" & `nshow' == `ncell' & `i' == `nshow' local lb "NR"
		if "`stat'" == "values" & `i' == `nshow' & `nshow' == `ncell' local lb "`lastlab'+"
		local xlab `"`xlab' `v' "`lb'""'
	}
	twoway (rcap ymin ymax xval, lcolor(gs10)) (rbar yp25 yp75 xval, barwidth(0.5) fcolor(gs14) lcolor(gs8)) ///
		(scatter ymedian xval, mcolor(black) msymbol(diamond)) (connected observed xval, lcolor(red) mcolor(red) msymbol(circle) lwidth(medthin)), ///
		xlabel(`xlab', angle(45) labsize(vsmall)) xtitle("`ttl'") ytitle("Proportion") title("`ttl'", size(small)) legend(off) nodraw name(`graphname', replace)
	restore
end

* estat ppcheck: posterior predictive check of the tie probabilities (expected values): observed values against the posterior mean predictions
capture program drop _nwlsm_ppcheck
program _nwlsm_ppcheck, rclass
	version 14
	syntax [, GRoups(integer 10) GRAPH NAME(string) NOList]
	if `groups' < 2 | `groups' > 50 {
		di as error "groups() must be between 2 and 50"
		exit 198
	}
	local fam "`e(family)'"
	if "`fam'" == "" local fam "bernoulli"
	local dir = ("`e(directed)'" == "directed")
	capture mata: st_numscalar("__nwlsm_Sck", rows(LsmSample("beta")))
	if _rc {
		di as error "the posterior draws of the last fit are no longer in memory; fit the model again"
		exit 498
	}
	scalar drop __nwlsm_Sck
	tempname T SM
	mata: LsmPPCheck(`groups', "`T'", "`SM'")
	matrix colnames `T' = dyads predicted observed
	local ng = rowsof(`T')
	if "`nolist'" == "" {
		di
		di as text "Posterior predictive check: observed " cond("`fam'" == "bernoulli", "ties", "values") " against the posterior mean " cond("`fam'" == "bernoulli", "tie probabilities", "expected values") " (`fam' family)"
		di as text "{hline 50}"
		di as text %8s "Group" %12s "Dyads" %15s "Predicted" %15s "Observed"
		di as text "{hline 50}"
		forvalues i = 1/`ng' {
			di as result %8.0f `i' %12.0f `T'[`i', 1] %15.4f `T'[`i', 2] %15.4f `T'[`i', 3]
		}
		di as text "{hline 50}"
		di as text "Observed dyads: " as result %8.0f `SM'[1, 1] _c
		if "`fam'" == "bernoulli" {
			di as text "  Brier score: " as result %7.4f `SM'[1, 2] as text "  AUC: " as result %7.4f `SM'[1, 3] as text "  mean log loss: " as result %7.4f `SM'[1, 4]
		}
		else {
			di as text "  RMSE: " as result %9.4f `SM'[1, 2] as text "  correlation: " as result %7.4f `SM'[1, 3] as text "  mean absolute error: " as result %9.4f `SM'[1, 4]
		}
		di as text "{p 0 0 2}Dyads are sorted by their posterior mean prediction and cut into `groups' groups of about equal size; a well calibrated model has the observed mean close to the predicted one in every group.{p_end}"
	}
	if "`graph'" != "" {
		if "`name'" == "" local name "ppcheck"
		preserve
		quietly drop _all
		quietly set obs `ng'
		mata: st_store(., st_addvar("double", "predicted"), st_matrix("`T'")[., 2])
		mata: st_store(., st_addvar("double", "observed"), st_matrix("`T'")[., 3])
		quietly summarize predicted, meanonly
		local lo = min(r(min), 0)
		quietly summarize observed, meanonly
		local lo = min(`lo', r(min))
		local hi = max(r(max), 1e-9)
		quietly summarize predicted, meanonly
		local hi = max(`hi', r(max))
		twoway (function y = x, range(`lo' `hi') lcolor(gs10) lpattern(dash)) (connected observed predicted, mcolor(navy) lcolor(navy) msymbol(circle)), ///
			xtitle("Posterior mean prediction (group mean)") ytitle("Observed (group mean)") title("Posterior predictive check", size(medium)) legend(off) name(`name', replace)
		restore
	}
	return matrix groups = `T'
	return scalar N = `SM'[1, 1]
	return scalar score = `SM'[1, 2]
	return scalar auc_or_corr = `SM'[1, 3]
end

* estat influence: the outcome equation of the latent space influence model, predictions and the actors who exert most influence
capture program drop _nwlsm_influence
program _nwlsm_influence, rclass
	version 14
	syntax [, PREDict(string) EXPosure(string) EXErted(string) TOP(integer 10) NOList GRaph NAME(string)]
	if "`e(influence)'" != "yes" {
		di as error "the model was fitted without influence()"
		exit 198
	}
	tempname T W PR EX
	matrix `T' = e(summary)
	local K : word count `e(covariates)'
	local q : word count `e(ocovariates)'
	local rb = `K' + `q'
	local n = e(N_nodes)
	capture matrix `PR' = e(pred)
	if _rc {
		di as error "the predictions were not stored (too many actors for a Stata matrix)"
		exit 498
	}
	matrix `W' = e(W)
	mata: __nwlsm_x = st_matrix("`PR'")
	mata: st_numscalar("__nwlsm_sst", sum((__nwlsm_x[., 1] :- mean(__nwlsm_x[., 1])) :^ 2))
	mata: st_numscalar("__nwlsm_sse", sum((__nwlsm_x[., 1] - __nwlsm_x[., 2]) :^ 2))
	mata: st_numscalar("__nwlsm_sse0", sum((__nwlsm_x[., 1] - (__nwlsm_x[., 2] - __nwlsm_x[., 4])) :^ 2))
	mata: st_matrix("`EX'", colsum(st_matrix("`W'"))')
	mata: st_numscalar("__nwlsm_mx", max(st_matrix("`EX'")))
	mata: st_numscalar("__nwlsm_ent", mean(1 :/ rowsum(st_matrix("`W'") :^ 2)))
	local r2 = 1 - scalar(__nwlsm_sse) / scalar(__nwlsm_sst)
	local r20 = 1 - scalar(__nwlsm_sse0) / scalar(__nwlsm_sst)
	local effn = scalar(__nwlsm_ent)
	scalar drop __nwlsm_sst __nwlsm_sse __nwlsm_sse0 __nwlsm_mx __nwlsm_ent
	local labs : rownames `W'
	if `"`predict'"' != "" | `"`exposure'"' != "" | `"`exerted'"' != "" {
		if _N != `n' {
			di as error "{p 0 0 2}predict(), exposure() and exerted() need a data set with one observation per actor in the order of the actors of {bf:`e(network)'}; the data set has `=_N' observations.{p_end}"
			exit 198
		}
		if `"`predict'"' != "" {
			confirm new variable `predict'
			quietly generate double `predict' = `PR'[_n, 2]
		}
		if `"`exposure'"' != "" {
			confirm new variable `exposure'
			quietly generate double `exposure' = `PR'[_n, 3]
		}
		if `"`exerted'"' != "" {
			confirm new variable `exerted'
			quietly generate double `exerted' = `EX'[_n, 1]
		}
	}
	if "`nolist'" == "" {
		di
		di as text "Latent space influence model: outcome " as result "`e(outcome)'" as text ", lagged outcome " as result "`e(lagvar)'" as text "; weights `e(weight)', " as result e(neighbors) as text " nearest actors"
		di as text "{hline 78}"
		di as text _col(25) %10s "Mean" %10s "Sd" %10s "MC s.e." %11s "Lower 95%" %10s "Upper 95%"
		di as text "{hline 78}"
		local k `K'
		foreach c in `e(ocovariates)' {
			local ++k
			di as text %-24s abbrev("`c'", 23) _col(25) as result %10.4f `T'[`k', 1] %10.4f `T'[`k', 2] %10.4f `T'[`k', 3] %11.4f `T'[`k', 4] %10.4f `T'[`k', 5]
		}
		tempname VT
		matrix `VT' = e(vsummary)
		local nv = rowsof(`VT')
		di as text %-24s "sigma2_y (error variance)" _col(25) as result %10.4f `VT'[`nv', 1] %10.4f `VT'[`nv', 2] %10.4f `VT'[`nv', 3] %11.4f `VT'[`nv', 4] %10.4f `VT'[`nv', 5]
		di as text "{hline 78}"
		di as text "Posterior probability that the influence effect is positive: " as result %6.4f e(prob_positive)
		di as text "R-squared of the posterior mean predictions: " as result %6.4f `r2' as text "; without the influence term: " as result %6.4f `r20' as text "; the influence term adds " as result %6.4f `r2' - `r20'
		di as text "Effective number of influencing actors per actor (inverse of the sum of squared posterior mean weights, mean): " as result %6.2f `effn'
		if `top' > 0 {
			mata: st_matrix("__nwlsm_ord", order(-st_matrix("`EX'"), 1))
			local m = min(`top', `n')
			di
			di as text "Actors who exert most influence (column sums of the posterior mean weights e(W))"
			di as text "{hline 62}"
			di as text %-22s "Actor" %20s "Influence exerted" %20s "Exposure received"
			di as text "{hline 62}"
			forvalues r = 1/`m' {
				local j = __nwlsm_ord[`r', 1]
				local lb : word `j' of `labs'
				if "`lb'" == "" local lb "`j'"
				di as text %-22s abbrev("`lb'", 21) as result %20.3f `EX'[`j', 1] %20.3f `PR'[`j', 3]
			}
			di as text "{hline 62}"
			capture matrix drop __nwlsm_ord
		}
	}
	if "`graph'" != "" {
		if "`name'" == "" local name "influence"
		preserve
		quietly {
			clear
			svmat double `PR', names(c)
			rename c1 obs
			rename c2 pred
		}
		twoway (scatter obs pred, mcolor(navy%70)) (function y = x, range(pred) lcolor(gs8) lpattern(dash)), ///
			xtitle("predicted outcome (posterior mean)") ytitle("observed outcome") title("Latent space influence model: fit") legend(off) ///
			name(`name', replace)
		restore
	}
	return scalar influence = `T'[`rb', 1]
	return scalar influence_sd = `T'[`rb', 2]
	return scalar influence_lower = `T'[`rb', 4]
	return scalar influence_upper = `T'[`rb', 5]
	return scalar prob_positive = e(prob_positive)
	return scalar r2 = `r2'
	return scalar r2_noinfluence = `r20'
	return scalar eff_neighbors = `effn'
	capture mata: mata drop __nwlsm_x
	return matrix exerted = `EX'
end

* estat bic: the BIC of the latent position model as latentnet's bic.ergmm (computed at the minimum Kullback-Leibler positions and effects)
capture program drop _nwlsm_bic
program _nwlsm_bic, rclass
	version 14
	syntax [, EFFobs(string) NOList]
	if "`e(family)'" != "bernoulli" & "`e(family)'" != "" {
		di as error "estat bic is available for binary networks (family bernoulli) only"
		exit 198
	}
	if "`e(latent)'" == "bilinear" {
		di as error "estat bic is not available for the bilinear model"
		exit 198
	}
	if "`e(influence)'" == "yes" {
		di as error "estat bic is not available for the influence model"
		exit 198
	}
	capture mata: st_numscalar("__nwlsm_Sck", rows(LsmSample("beta")))
	if _rc {
		di as error "the posterior draws of the last fit are no longer in memory; fit the model again"
		exit 498
	}
	scalar drop __nwlsm_Sck
	local eff 0
	local effname "ties"
	if "`effobs'" != "" {
		local effname = lower("`effobs'")
		local n = e(N_nodes)
		local dir = ("`e(directed)'" == "directed")
		if "`effname'" == "ties" local eff 0
		else if "`effname'" == "dyads" local eff = cond(`dir', `n' * (`n' - 1), `n' * (`n' - 1) / 2)
		else if "`effname'" == "actors" local eff `n'
		else {
			di as error "effobs() must be ties, dyads or actors"
			exit 198
		}
	}
	tempname E R
	local n = e(N_nodes)
	capture matrix `E' = e(mkl_effects)
	if _rc {
		capture confirm matrix e(Z)
		if _rc {
			di as error "the minimum Kullback-Leibler estimates were not stored (too many actors for a Stata matrix)"
			exit 498
		}
		matrix `E' = J(`n', 3, 0)
	}
	local G = e(clusters)
	local effl `e(effects)'
	mata: st_matrix("`R'", LsmBicFit(st_matrix("`E'"), `G', `eff'))
	matrix colnames `R' = Y Z sender receiver sociality overall p loglik Gused effobs
	local ob = `R'[1, 6]
	if "`nolist'" == "" {
		di
		di as text "BIC of the latent position model (as latentnet's bic.ergmm, at the minimum Kullback-Leibler estimates)"
		di as text "{hline 44}"
		di as text %-24s "Component" %20s "BIC"
		di as text "{hline 44}"
		di as text %-24s "Ties given positions (Y)" as result %20.3f `R'[1, 1]
		di as text %-24s "Positions (Z)" as result %20.3f `R'[1, 2]
		foreach t in sender receiver sociality {
			local k = cond("`t'" == "sender", 3, cond("`t'" == "receiver", 4, 5))
			if `: list t in effl' di as text %-24s "Random effect: `t'" as result %20.3f `R'[1, `k']
		}
		di as text "{hline 44}"
		di as text %-24s "Overall" as result %20.3f `ob'
		di as text "{hline 44}"
		di as text "Latent clusters: " as result `G' as text "; coefficients: " as result `R'[1, 7] as text "; effective number of observations (" "`effname'" "): " as result `R'[1, 10]
		di as text "{p 0 0 2}Smaller is better. The BIC compares models fitted to the same network; it is not certain that it is appropriate for choosing the dimension, the random effects or between a clustered and an unclustered model.{p_end}"
	}
	return scalar bic = `ob'
	return scalar bic_y = `R'[1, 1]
	return scalar bic_z = `R'[1, 2]
	return scalar bic_sender = `R'[1, 3]
	return scalar bic_receiver = `R'[1, 4]
	return scalar bic_sociality = `R'[1, 5]
	return scalar p = `R'[1, 7]
	return scalar loglik_mle = `R'[1, 8]
	return scalar effobs = `R'[1, 10]
	return scalar clusters_used = `R'[1, 9]
	return matrix components = `R'
end
