* Agreement with the R package latentnet (reference results in dev/ref, made by dev/ref_latentnet.R): Sampson monks (directed, d = 2),
* a random undirected network with a node covariate (absdiff), n = 40, d = 2 and d = 3. Four seeds per model. Criteria:
*  - posterior mean of every coefficient (average over the seeds) within 3.5 Monte Carlo standard errors of the latentnet value, where the
*    standard error combines the batch-means standard error of our average (4 chains) and one of the same size for the latentnet chain;
*  - posterior mean log likelihood within 2 posterior standard deviations of the latentnet value;
*  - correlation of the upper triangles of the posterior mean distance matrices at least 0.98;
*  - correlation of the pairwise distances of the minimum-KL positions with those of latentnet at least 0.90 (the objective of the minimum-KL fit has several local optima, so that the
*    two programs can end in different ones in a sparse network; nwlsm starts from the best 5 draws and keeps the best optimum).
do cscripts/nwlsm_testlib.do
global NS 4
capture program drop lsmt_cmp
program lsmt_cmp
	syntax name(name=case), N(integer) DIR(integer) D(integer) [X]
	local ref "$LSM_REPO/dev/ref"
	lsmt_mknet net_`case', edges("`ref'/`case'_edges.csv") n(`n') dir(`dir')
	tempname RB RL MD
	* reference numbers
	preserve
	import delimited "`ref'/`case'_beta.csv", clear varnames(1)
	mkmat mean sd, matrix(`RB')
	import delimited "`ref'/`case'_llk.csv", clear varnames(1)
	local rllk = llk_mean[1]
	local rllksd = llk_sd[1]
	import delimited "`ref'/`case'_meandist.csv", clear varnames(1)
	mkmat *, matrix(`MD')
	import delimited "`ref'/`case'_mklZ.csv", clear varnames(1)
	mkmat *, matrix(__RZ)
	restore
	local opt
	if "`x'" != "" {
		preserve
		import delimited "`ref'/random40_x.csv", clear varnames(1)
		mkmat x, matrix(__x)
		mata: __x = st_matrix("__x"); st_matrix("__AD", abs(__x * J(1, `n', 1) - J(`n', 1, 1) * __x'))
		restore
		local opt "edgecov(__AD)"
	}
	local K = rowsof(`RB')
	matrix SB = J($NS, `K', .)
	matrix SM = J($NS, `K', .)
	local llsum 0
	tempname DS
	matrix `DS' = J(`n', `n', 0)
	forvalues s = 1/$NS {
		quietly nwlsm net_`case', dim(`d') `opt' seed(`=1000 + `s'')
		matrix b = e(b)
		matrix T = e(summary)
		forvalues k = 1/`K' {
			matrix SB[`s', `k'] = b[1, `k']
			matrix SM[`s', `k'] = T[`k', 3]
		}
		if `s' == 1 matrix MKLZ = e(Z)
		local llsum = `llsum' + e(ll_mean)
		matrix `DS' = `DS' + e(D) / $NS
	}
	local llavg = `llsum' / $NS
	local ok 1
	forvalues k = 1/`K' {
		mata: st_local("avg", strofreal(mean(st_matrix("SB")[., `k'])))
		mata: st_local("se", strofreal(sqrt(mean(st_matrix("SM")[., `k'] :^ 2) / $NS + mean(st_matrix("SM")[., `k'] :^ 2))))
		local diff = `avg' - `RB'[`k', 1]
		local z = abs(`diff') / `se'
		di as text "case `case' coef `k': nwlsm " %8.4f `avg' "  latentnet " %8.4f `RB'[`k', 1] "  (posterior sd " %6.4f `RB'[`k', 2] ")  diff/MCse = " %5.2f `z'
		if `z' > 3.5 local ok 0
	}
	local lz = abs(`llavg' - `rllk') / `rllksd'
	di as text "case `case' log likelihood: nwlsm " %9.3f `llavg' "  latentnet " %9.3f `rllk' "  diff/sd = " %5.2f `lz'
	if `lz' > 2 local ok 0
	mata: st_local("cor", strofreal(lsmt_cor_ut(st_matrix("`DS'"), st_matrix("`MD'")), "%9.5f"))
	di as text "case `case' correlation of posterior mean distances: " `cor'
	if `cor' < 0.98 local ok 0
	mata: st_local("corz", strofreal(lsmt_cor_ut(lsmt_pd(st_matrix("MKLZ")), lsmt_pd(st_matrix("__RZ"))), "%9.5f"))
	di as text "case `case' correlation of the pairwise distances of the minimum-KL positions with latentnet's: " `corz'
	if `corz' < 0.90 local ok 0
	if !`ok' {
		di as error "FAIL case `case'"
		exit 9
	}
	di as result "PASS case `case'"
end
mata:
real scalar lsmt_cor_ut(real matrix A, real matrix B)
{
	real colvector a, b
	real scalar n, i, j, k

	n = rows(A)
	a = J(n * (n - 1) / 2, 1, .); b = a; k = 0
	for (i = 1; i < n; i++) for (j = i + 1; j <= n; j++) {
		k++; a[k] = A[i, j]; b[k] = B[i, j]
	}
	return(correlation((a, b))[1, 2])
}
end
lsmt_cmp sampson_d2, n(18) dir(1) d(2)
lsmt_cmp random40_d2, n(40) dir(0) d(2) x
lsmt_cmp random40_d3, n(40) dir(0) d(3) x
display "PASS test_nwlsm_latentnet: Sampson d=2, random40 d=2 and d=3 agree with latentnet"
