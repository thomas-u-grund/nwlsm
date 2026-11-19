* Agreement with the R package latentnet for the models of version 0.3.0: sender and receiver effects (directed), sociality (undirected and
* directed), a node covariate, the bilinear model. Reference results in dev/ref (two latentnet chains per case, made by dev/ref_latentnet.R); six
* nwlsm chains per case (seeds 1001 to 1006). Criteria:
*  - posterior mean of every coefficient and variance parameter (average over the seeds) within 4 standard errors of the reference mean (average of
*    the two latentnet chains), where the standard error combines the batch-means error of our average, a standard error for the reference taken
*    from the difference of its two chains (at least 1/7 of the posterior sd), 
*  - posterior mean log likelihood within 2 posterior sd of the reference;
*  - correlation of the posterior mean distances (Euclidean models) with the reference at least 0.98;
*  - correlation of the posterior mean tie probabilities (off-diagonal) with the reference at least 0.995 and mean absolute difference below 0.01;
*  - correlation of the posterior mean effects with the reference at least 0.97;
*  - correlation of the pairwise distances of the minimum-KL positions (Euclidean models) with those of latentnet at least 0.98.
do cscripts/nwlsm_testlib.do
global NS 6
capture program drop lsmt_cmp2
program lsmt_cmp2
	syntax name(name=case), N(integer) DIR(integer) D(integer) [EFFects(string) X(string) BILinear]
	local ref "$LSM_REPO/dev/ref"
	lsmt_mknet net_`case', edges("`ref'/`case'_edges.csv") n(`n') dir(`dir')
	tempname RB RB2 RV RV2 DS PS ES
	preserve
	import delimited "`ref'/`case'_beta.csv", clear varnames(1)
	mkmat mean sd, matrix(`RB')
	import delimited "`ref'/`case'_c2_beta.csv", clear varnames(1)
	mkmat mean sd, matrix(`RB2')
	import delimited "`ref'/`case'_vars.csv", clear varnames(1)
	mkmat mean sd, matrix(`RV')
	import delimited "`ref'/`case'_c2_vars.csv", clear varnames(1)
	mkmat mean sd, matrix(`RV2')
	import delimited "`ref'/`case'_llk.csv", clear varnames(1)
	local rllk = llk_mean[1]
	local rllksd = llk_sd[1]
	import delimited "`ref'/`case'_meandist.csv", clear varnames(1)
	mkmat *, matrix(__MD1)
	import delimited "`ref'/`case'_c2_meandist.csv", clear varnames(1)
	mkmat *, matrix(__MD2)
	import delimited "`ref'/`case'_P.csv", clear varnames(1)
	mkmat *, matrix(__RP)
	import delimited "`ref'/`case'_mklZ.csv", clear varnames(1)
	mkmat *, matrix(__RZ)
	if "`effects'" != "" {
		import delimited "`ref'/`case'_effects.csv", clear varnames(1)
		mkmat *, matrix(__RE)
	}
	restore
	local opt "`effects' `bilinear'"
	if "`x'" != "" {
		preserve
		import delimited "`ref'/`x'", clear varnames(1)
		mkmat x, matrix(__x)
		mata: __x = st_matrix("__x"); st_matrix("__AD", abs(__x * J(1, `n', 1) - J(`n', 1, 1) * __x'))
		restore
		local opt "`opt' edgecov(__AD)"
	}
	local K = rowsof(`RB')
	local V = rowsof(`RV')
	matrix SB = J($NS, `K', .)
	matrix SBm = J($NS, `K', .)
	matrix SV = J($NS, `V', .)
	matrix SVm = J($NS, `V', .)
	local llsum 0
	matrix `DS' = J(`n', `n', 0)
	matrix `PS' = J(`n', `n', 0)
	local ne : word count `effects'
	if `ne' > 0 matrix `ES' = J(`n', `ne', 0)
	forvalues s = 1/$NS {
		quietly nwlsm net_`case', dim(`d') `opt' seed(`=1000 + `s'')
		matrix b = e(b)
		matrix T = e(summary)
		matrix VT = e(vsummary)
		forvalues k = 1/`K' {
			matrix SB[`s', `k'] = b[1, `k']
			matrix SBm[`s', `k'] = T[`k', 3]
		}
		forvalues k = 1/`V' {
			matrix SV[`s', `k'] = VT[`k', 1]
			matrix SVm[`s', `k'] = VT[`k', 3]
		}
		local llsum = `llsum' + e(ll_mean)
		matrix `DS' = `DS' + e(D) / $NS
		matrix `PS' = `PS' + e(P) / $NS
		local j 0
		foreach ef of local effects {
			local ++j
			matrix __e = e(`ef')
			forvalues i = 1/`n' {
				matrix `ES'[`i', `j'] = `ES'[`i', `j'] + __e[`i', 1] / $NS
			}
		}
		if `s' == 1 matrix MKLZ = e(Z)
	}
	local ok 1
	forvalues k = 1/`K' {
		mata: st_local("avg", strofreal(mean(st_matrix("SB")[., `k'])))
		mata: st_local("se", strofreal(sqrt(mean(st_matrix("SBm")[., `k'] :^ 2) / $NS)))
		local ref = (`RB'[`k', 1] + `RB2'[`k', 1]) / 2
		local seref = max(abs(`RB'[`k', 1] - `RB2'[`k', 1]) / 2, `RB'[`k', 2] / 7)
		local z = abs(`avg' - `ref') / sqrt(`se'^2 + `seref'^2)
		di as text "case `case' coef `k': nwlsm " %8.4f `avg' "  latentnet " %8.4f `ref' "  (posterior sd " %6.4f `RB'[`k', 2] ")  diff/se = " %5.2f `z'
		if `z' > 4 local ok 0
	}
	forvalues k = 1/`V' {
		mata: st_local("avg", strofreal(mean(st_matrix("SV")[., `k'])))
		mata: st_local("se", strofreal(sqrt(mean(st_matrix("SVm")[., `k'] :^ 2) / $NS)))
		local ref = (`RV'[`k', 1] + `RV2'[`k', 1]) / 2
		local seref = max(abs(`RV'[`k', 1] - `RV2'[`k', 1]) / 2, `RV'[`k', 2] / 7)
		local z = abs(`avg' - `ref') / sqrt(`se'^2 + `seref'^2)
		local nmv : word `k' of `: rownames VT'
		di as text "case `case' variance `nmv': nwlsm " %8.4f `avg' "  latentnet " %8.4f `ref' "  (posterior sd " %6.4f `RV'[`k', 2] ")  diff/se = " %5.2f `z'
		if `z' > 4 local ok 0
	}
	local llavg = `llsum' / $NS
	local lz = abs(`llavg' - `rllk') / `rllksd'
	di as text "case `case' log likelihood: nwlsm " %9.3f `llavg' "  latentnet " %9.3f `rllk' "  diff/sd = " %5.2f `lz'
	if `lz' > 2 local ok 0
	if "`bilinear'" == "" {
		mata: st_local("cor", strofreal(lsmt_cor_ut(st_matrix("`DS'"), (st_matrix("__MD1") + st_matrix("__MD2")) / 2), "%9.5f"))
		di as text "case `case' correlation of posterior mean distances: " `cor'
		if `cor' < 0.98 local ok 0
	}
	* the minimum-KL estimate: its objective (likelihood of the posterior mean tie probabilities) equals the objective recomputed from the
	* returned positions, coefficients and effects, and is not below the objective at the posterior means
	mata: st_local("o1", strofreal(lsmt_mklobj(st_matrix("e(P)"), `dir', st_matrix("e(mkl_b)"), st_matrix("e(Z)"), st_matrix("e(mkl_effects)"), "`effects'", "`bilinear'", `n', 0), "%20.10f"))
	local o0 = e(mkl_obj)
	di as text "case `case' MKL objective: plugin " %12.6f `o0' "  recomputed " %12.6f `o1'
	if abs(`o0' - `o1') > 1e-6 local ok 0
	mata: st_local("o2", strofreal(lsmt_mklobj(st_matrix("e(P)"), `dir', st_matrix("e(b)"), st_matrix("e(Zmean)"), lsmt_effmeans("`effects'", `n'), "`effects'", "`bilinear'", `n', 1), "%20.10f"))
	di as text "case `case' objective at the posterior means (positions aligned, not optimised): " %12.6f `o2'
	if `o0' < `o2' - 1e-6 local ok 0
	mata: st_local("corp", strofreal(lsmt_cor_ut(st_matrix("`PS'"), st_matrix("__RP")), "%9.5f"))
	mata: st_local("mad", strofreal((sum(abs(st_matrix("`PS'") - st_matrix("__RP"))) - sum(abs(diagonal(st_matrix("`PS'") - st_matrix("__RP"))))) / (`n' * (`n' - 1)), "%9.5f"))
	di as text "case `case' posterior mean tie probabilities: correlation " `corp' ", mean absolute difference " `mad'
	if `corp' < 0.995 | `mad' > 0.01 local ok 0
	if `ne' > 0 {
		forvalues j = 1/`ne' {
			mata: st_local("cr", strofreal(correlation((st_matrix("`ES'")[., `j'], st_matrix("__RE")[., `j']))[1, 2], "%9.5f"))
			local ef : word `j' of `effects'
			di as text "case `case' correlation of the posterior mean `ef' effects: " `cr'
			if `cr' < 0.97 local ok 0
		}
	}
	if !`ok' {
		di as error "FAIL case `case'"
		exit 9
	}
	di as result "PASS case `case'"
end
mata:
mata set matastrict off
// effects of the last fit: posterior means (n x 3 with zero columns for the absent ones)
real matrix lsmt_effmeans(string scalar effects, real scalar n)
{
	real matrix E, M
	string rowvector nm
	real scalar j

	E = J(n, 3, 0)
	nm = tokens(effects)
	for (j = 1; j <= cols(nm); j++) {
		M = st_matrix("e(" + nm[j] + ")")
		E[., (nm[j] == "sender") + 2 * (nm[j] == "receiver") + 3 * (nm[j] == "sociality")] = M[., 1]
	}
	return(E)
}
// the MKL objective sum over dyads of EY eta - log(1 + exp(eta)) at given parameters (the effects: n x 3 columns sender, receiver, sociality;
// mkl_effects has three columns, effmeans too)
real scalar lsmt_mklobj(real matrix P, real scalar dir, real rowvector b, real matrix Z, real matrix E, string scalar effects, string scalar bil, real scalar n, real scalar ispost)
{
	real scalar i, j, eta, f, lt, hs, hr, ho
	hs = strpos(" " + effects + " ", " sender ") > 0
	hr = strpos(" " + effects + " ", " receiver ") > 0
	ho = strpos(" " + effects + " ", " sociality ") > 0
	f = 0
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j | (!dir & j < i)) continue
		eta = b[1]
		if (cols(b) > 1) eta = eta + b[2] * lsmt_cov(i, j)
		if (hs) eta = eta + E[i, 1]
		if (hr) eta = eta + E[j, 2]
		if (ho) eta = eta + E[i, 3] + E[j, 3]
		if (bil != "") lt = -sum(Z[i, .] :* Z[j, .])
		else lt = sqrt(sum((Z[i, .] - Z[j, .]) :^ 2))
		eta = eta - lt
		f = f + P[i, j] * eta - ln(1 + exp(eta))
	}
	return(f)
}
real scalar lsmt_cov(real scalar i, real scalar j)
{
	return(st_matrix("__AD")[i, j])
}
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
lsmt_cmp2 sampson_sr, n(18) dir(1) d(2) effects(sender receiver)
lsmt_cmp2 sampson_soc, n(18) dir(1) d(2) effects(sociality)
lsmt_cmp2 soc40, n(40) dir(0) d(2) effects(sociality) x(soc40_x.csv)
lsmt_cmp2 sr40, n(40) dir(1) d(2) effects(sender receiver) x(sr40_x.csv)
lsmt_cmp2 bil40, n(40) dir(0) d(2) bilinear
display "PASS test_nwlsm_effects: sender, receiver and sociality effects, the MKL positions and the bilinear model agree with latentnet"
