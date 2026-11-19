* Influence model, recovery and coverage: data simulated from the exact model (ties by latent distance, outcome by the weighted lagged outcomes of the
* G nearest actors) with known parameters; REPS replications per scenario (default 30; env NWLSM_REPS), three scenarios (the wrappers test_nwlsm_infl_rec_n50, _n100, _all run one each, so that they can run in parallel). For each replication the 95% posterior
* interval of the network intercept, the outcome coefficients (intercept, lagged outcome, covariate, influence) and the error variance is compared with the truth;
* the posterior mean distances are compared with the true distances (correlation) and the positions with the true positions up to rotation, reflection and translation.
do cscripts/nwlsm_testlib.do
do cscripts/nwlsm_infl_sim.do
local reps : env NWLSM_REPS
if "`reps'" == "" local reps 30
local only "$LSMI_SCEN"
capture program drop lsmi_scen
program lsmi_scen
	syntax, N(integer) G(integer) WT(integer) REPS(integer) SEED(integer) NAME(string) [DIR(integer 1) B2(real 1)]
	local wtn = cond(`wt' == 0, "invdist", "expdist")
	tempname R
	matrix `R' = J(`reps', 9, .)
	forvalues r = 1/`reps' {
		lsmi_sim, n(`n') dir(`dir') d(2) a0(1.5) b0(0.5) b1(0.6) b2(`b2') b3(-0.4) s2(0.25) g(`g') wt(`wt') seed(`=`seed' + `r'') zsd(2.5)
		nwset __a*, name(sc) `=cond(`dir', "directed", "undirected")' replace
		quietly nwlsm sc, dim(2) influence(y ylag x) neighbors(`g') weight(`wtn') seed(`=`seed' + 1000 + `r'') nsample(2000) burnin(3000) interval(5) threads(1)
		matrix T = e(summary)
		matrix V = e(vsummary)
		local nv = rowsof(V)
		local truth 1.5 0.5 0.6 -0.4 `b2'
		* rows of e(summary): network _cons, outcome _cons, ylag, x, influence
		local cov 0
		matrix `R'[`r', 1] = (1.5 >= T[1, 4] & 1.5 <= T[1, 5])
		matrix `R'[`r', 2] = (0.5 >= T[2, 4] & 0.5 <= T[2, 5])
		matrix `R'[`r', 3] = (0.6 >= T[3, 4] & 0.6 <= T[3, 5])
		matrix `R'[`r', 4] = (-0.4 >= T[4, 4] & -0.4 <= T[4, 5])
		matrix `R'[`r', 5] = (`b2' >= T[5, 4] & `b2' <= T[5, 5])
		matrix `R'[`r', 6] = (0.25 >= V[`nv', 4] & 0.25 <= V[`nv', 5])
		matrix `R'[`r', 7] = T[5, 1] - `b2'
		mata: st_matrix("__c", (lsmi_cor_ut(st_matrix("e(D)"), __lsmi_D), lsmi_procrmse(st_matrix("e(Z)"), __lsmi_Z)))
		matrix `R'[`r', 8] = __c[1, 1]
		matrix `R'[`r', 9] = __c[1, 2]
	}
	mata: st_matrix("__m", mean(st_matrix("`R'")))
	mata: st_matrix("__sd", sqrt(diagonal(variance(st_matrix("`R'"))))')
	di as text "SCENARIO `name': n=`n' G=`g' weight=`wtn' dir=`dir' b2=`b2', `reps' replications"
	di as text "  coverage of 95% intervals: network intercept " %5.3f __m[1,1] "  outcome intercept " %5.3f __m[1,2] "  lag " %5.3f __m[1,3] "  covariate " %5.3f __m[1,4] "  influence " %5.3f __m[1,5] "  error variance " %5.3f __m[1,6]
	di as text "  influence: mean bias " %7.4f __m[1,7] "  sd of the posterior mean over replications " %7.4f __sd[1,7] "  (mean posterior mean " %6.3f `b2' + __m[1,7] ")"
	di as text "  correlation of true and posterior mean distances " %6.3f __m[1,8] "  Procrustes rmse of the positions " %6.3f __m[1,9]
	matrix `R' = `R'
	local tot = 0
	forvalues j = 1/6 {
		local tot = `tot' + __m[1, `j']
	}
	* coverage of at least 80% for every parameter (30 replications: the 95% interval misses 6 of 30 with probability about 0.0001)
	forvalues j = 1/6 {
		if __m[1, `j'] < 0.80 {
			di as error "FAIL scenario `name': coverage of parameter `j' is " %5.3f __m[1, `j']
			exit 9
		}
	}
	if abs(__m[1, 7]) > 0.35 {
		di as error "FAIL scenario `name': bias of the influence effect"
		exit 9
	}
	di as result "PASS scenario `name'"
end
if "`only'" == "" | "`only'" == "1" lsmi_scen, n(50) g(10) wt(0) reps(`reps') seed(5000) name(n50_G10_invdist)
if "`only'" == "" | "`only'" == "2" lsmi_scen, n(100) g(10) wt(0) reps(`reps') seed(6000) name(n100_G10_invdist)
if "`only'" == "" | "`only'" == "3" lsmi_scen, n(50) g(0) wt(1) reps(`reps') seed(7000) name(n50_all_expdist) dir(0)
display "PASS nwlsm_infl_recovery"
