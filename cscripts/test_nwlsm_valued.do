* Valued networks (nwlsm 0.5.0): agreement with the R package latentnet for the families Poisson, binomial (trials) and normal.
* Reference results in dev/ref (two latentnet chains per case, made by dev/ref_latentnet.R); six nwlsm chains per case (seeds 1001 to 1006). Criteria:
*  - posterior mean of every coefficient and variance parameter (also the dispersion) within 4 standard errors of the reference mean (average of the two
*    latentnet chains), the standard error combining the batch-means error of our average and one for the reference taken from the difference of its chains
*    (at least 1/7 of the posterior sd);
*  - posterior mean log likelihood within 2 posterior sd of the reference;
*  - correlation of the posterior mean distances with the reference at least 0.98;
*  - correlation of the posterior mean expected values with the reference at least 0.995 and mean absolute difference below 2% of their mean;
*  - correlation of the posterior mean effects with the reference at least 0.97;
*  - the minimum-KL objective equals its recomputation from the returned parameters, is not below the objective at latentnet's minimum-KL
*    parameters evaluated on the same expected values (the optimum of ours is at least as good), and not below the objective at the posterior means.
do cscripts/nwlsm_testlib.do
set linesize 255
global NS 6
capture program drop lsmt_cmpv
program lsmt_cmpv
	syntax name(name=case), N(integer) DIR(integer) D(integer) FAM(string) [EFFects(string) X(string) TRIals(string) DISPVar(string) DISPDf(string) NOMKLcmp]
	local ref "$LSM_REPO/dev/ref"
	lsmt_mkvnet net_`case', y("`ref'/`case'_y.csv") dir(`dir')
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
	import delimited "`ref'/`case'_mklpar.csv", clear varnames(1)
	quietly destring *, replace force
	mkmat *, matrix(__RM)
	if "`effects'" != "" {
		import delimited "`ref'/`case'_effects.csv", clear varnames(1)
		mkmat *, matrix(__RE)
	}
	restore
	* the V rows of the reference variance file are Z.var, the effect variances and the dispersion; ours: sigma2, effects, dispersion (same order)
	local opt "family(`fam') `effects'"
	if "`x'" != "" {
		preserve
		import delimited "`ref'/`x'", clear varnames(1)
		mkmat x, matrix(__x)
		mata: __x = st_matrix("__x"); st_matrix("__AD", abs(__x * J(1, `n', 1) - J(`n', 1, 1) * __x'))
		restore
		local opt "`opt' edgecov(__AD)"
	}
	else matrix __AD = J(`n', `n', 0)
	if "`trials'" != "" {
		if "`trials'" == "mat" {
			preserve
			import delimited "`ref'/`case'_tr.csv", clear varnames(1)
			mkmat *, matrix(__TR)
			restore
			local opt "`opt' trials(__TR)"
		}
		else local opt "`opt' trials(`trials')"
	}
	if "`dispvar'" != "" local opt "`opt' dispvar(`dispvar') dispdf(`dispdf')"
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
		if `s' == 1 {
			matrix MKLZ = e(Z)
			matrix MKLB = e(mkl_b)
			matrix MKLE = e(mkl_effects)
			matrix MKLV = e(vsummary)
			local o0 = e(mkl_obj)
			matrix ZMEAN = e(Zmean)
			matrix BMEAN = e(b)
			matrix PMKL = e(P)
		}
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
	mata: st_local("cor", strofreal(lsmt_cor_ut(st_matrix("`DS'"), (st_matrix("__MD1") + st_matrix("__MD2")) / 2), "%9.5f"))
	di as text "case `case' correlation of posterior mean distances: " `cor'
	if `cor' < 0.98 local ok 0
	mata: st_local("corp", strofreal(lsmt_cor_ut(st_matrix("`PS'"), st_matrix("__RP")), "%9.5f"))
	mata: st_local("mad", strofreal(lsmt_mad(st_matrix("`PS'"), st_matrix("__RP")), "%9.5f"))
	di as text "case `case' posterior mean expected values: correlation " `corp' ", mean absolute difference relative to the mean " `mad'
	if `corp' < 0.995 | `mad' > 0.02 local ok 0
	if `ne' > 0 {
		forvalues j = 1/`ne' {
			mata: st_local("cr", strofreal(correlation((st_matrix("`ES'")[., `j'], st_matrix("__RE")[., `j']))[1, 2], "%9.5f"))
			local ef : word `j' of `effects'
			di as text "case `case' correlation of the posterior mean `ef' effects: " `cr'
			if `cr' < 0.97 local ok 0
		}
	}
	* minimum-KL objective on the expected values PMKL of the first chain
	local dm = 1
	if "`fam'" == "normal" local dm = MKLV[rowsof(MKLV), 1]
	mata: lsmt_fam = "`fam'"
	mata: st_local("o1", strofreal(lsmt_mklobjv(st_matrix("PMKL"), `dir', st_matrix("MKLB"), st_matrix("MKLZ"), st_matrix("MKLE"), "`effects'", `n', `dm', "`trials'"), "%20.10f"))
	di as text "case `case' MKL objective: plugin " %14.6f `o0' "  recomputed " %14.6f `o1'
	if abs(`o0' - `o1') > 1e-6 * (1 + abs(`o0')) local ok 0
	* latentnet's minimum-KL parameters on the same expected values
	mata: st_local("o3", strofreal(lsmt_mklobjv_ref(st_matrix("PMKL"), `dir', st_matrix("__RM"), st_matrix("__RZ"), "`effects'", `n', `dm', "`trials'"), "%20.10f"))
	di as text "case `case' MKL objective at latentnet's minimum-KL parameters (same expected values): " %14.6f `o3' "  (nwlsm: " %14.6f `o0' ")"
	if `o0' < `o3' - 1e-6 * (1 + abs(`o3')) local ok 0
	mata: st_local("o2", strofreal(lsmt_mklobjv(st_matrix("PMKL"), `dir', st_matrix("BMEAN"), st_matrix("ZMEAN"), lsmt_effmeans("`effects'", `n'), "`effects'", `n', `dm', "`trials'", 1), "%20.10f"))
	di as text "case `case' objective at the posterior means (positions aligned, not optimised): " %14.6f `o2'
	if `o0' < `o2' - 1e-6 * (1 + abs(`o2')) local ok 0
	if !`ok' {
		di as error "FAIL case `case'"
		exit 9
	}
	di as result "PASS case `case'"
end
mata:
mata set matastrict off
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
// objective -sum klterm over the observed dyads at given parameters; fam in lsmt_fam; trials: "" none, "mat" the matrix __TR, else a number
real scalar lsmt_mklobjv_core(real matrix P, real scalar dir, real rowvector b, real matrix Z, real matrix E, string scalar effects, real scalar n, real scalar dm, string scalar trials)
{
	real scalar i, j, eta, f, lt, hs, hr, ho, tr
	external lsmt_fam
	hs = strpos(" " + effects + " ", " sender ") > 0
	hr = strpos(" " + effects + " ", " receiver ") > 0
	ho = strpos(" " + effects + " ", " sociality ") > 0
	f = 0
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j | (!dir & j < i)) continue
		eta = b[1]
		if (cols(b) > 1) eta = eta + b[2] * st_matrix("__AD")[i, j]
		if (hs) eta = eta + E[i, 1]
		if (hr) eta = eta + E[j, 2]
		if (ho) eta = eta + E[i, 3] + E[j, 3]
		lt = sqrt(sum((Z[i, .] - Z[j, .]) :^ 2))
		eta = eta - lt
		if (trials == "") tr = 1
		else if (trials == "mat") tr = st_matrix("__TR")[i, j]
		else tr = strtoreal(trials)
		if (lsmt_fam == "poisson") f = f - (exp(eta) - P[i, j] * eta)
		else if (lsmt_fam == "binomial") f = f - (tr * ln(1 + exp(eta)) - P[i, j] * eta)
		else f = f - (0.5 * eta * eta - P[i, j] * eta) / dm
	}
	return(f)
}
real scalar lsmt_mklobjv(real matrix P, real scalar dir, real rowvector b, real matrix Z, real matrix E, string scalar effects, real scalar n, real scalar dm, string scalar trials, | real scalar ispost)
{
	return(lsmt_mklobjv_core(P, dir, b, Z, E, effects, n, dm, trials))
}
// latentnet's parameters: RM columns beta, then the effects it fitted (sender, receiver, sociality; NA padded)
real scalar lsmt_mklobjv_ref(real matrix P, real scalar dir, real matrix RM, real matrix Z, string scalar effects, real scalar n, real scalar dm, string scalar trials)
{
	real rowvector b
	real matrix E
	real scalar j, c
	string rowvector nm

	b = select(RM[., 1], RM[., 1] :< .)'
	E = J(n, 3, 0)
	nm = tokens(effects)
	c = 1
	for (j = 1; j <= cols(nm); j++) {
		c++
		E[., (nm[j] == "sender") + 2 * (nm[j] == "receiver") + 3 * (nm[j] == "sociality")] = RM[1..n, c]
	}
	return(lsmt_mklobjv_core(P, dir, b, Z, E, effects, n, dm, trials))
}
real scalar lsmt_mad(real matrix A, real matrix B)
{
	real scalar n

	n = rows(A)
	return((sum(abs(A - B)) - sum(abs(diagonal(A - B)))) / (n * (n - 1)) / (sum(abs(B)) / (n * n)))
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
lsmt_cmpv pois30, n(30) dir(1) d(2) fam(poisson) x(pois30_x.csv)
lsmt_cmpv pois30s, n(30) dir(0) d(2) fam(poisson) effects(sociality)
lsmt_cmpv bin30, n(30) dir(1) d(2) fam(binomial) effects(sender receiver) x(bin30_x.csv) trials(6)
lsmt_cmpv bin25m, n(25) dir(0) d(2) fam(binomial) trials(mat)
lsmt_cmpv gau30, n(30) dir(1) d(2) fam(normal) x(gau30_x.csv) dispvar(1) dispdf(3)
lsmt_cmpv gau30s, n(30) dir(0) d(2) fam(normal) effects(sociality) dispvar(0.5) dispdf(4)
display "PASS test_nwlsm_valued: Poisson, binomial and normal families agree with latentnet"
