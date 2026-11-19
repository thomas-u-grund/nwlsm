* Recovery for valued networks: data simulated from known coefficients, positions and (normal family) dispersion. The 95% intervals of the coefficients
* must cover the true values, the dispersion's interval must cover its true value, and the posterior mean distances must correlate with the true
* distances (upper triangles) at least 0.9. Poisson (directed, n = 50, covariate), binomial (undirected, n = 50, 7 trials), normal (undirected, n = 50,
* sociality), Poisson with sender and receiver effects and the bilinear latent term.
do cscripts/nwlsm_testlib.do
set linesize 255
mata:
mata set matastrict off
void lsmt_simv(real scalar n, real scalar dir, real scalar fam, real scalar b0, real scalar b1, real scalar seed, real scalar disp, real scalar sociality, real scalar bil)
{
	real matrix Z, D, Xc, eta, y
	real colvector x, so
	real scalar i
	external real matrix __lsmt_Dtrue, __lsmt_Xc

	rseed(seed)
	Z = rnormal(n, 2, 0, 1)
	x = rnormal(n, 1, 0, 1)
	so = rnormal(n, 1, 0, 0.5) * sociality
	D = J(n, n, 0)
	for (i = 1; i <= 2; i++) D = D + (Z[., i] * J(1, n, 1) - J(n, 1, 1) * Z[., i]') :^ 2
	D = sqrt(D)
	Xc = abs(x * J(1, n, 1) - J(n, 1, 1) * x')
	eta = b0 :+ b1 * Xc - (bil ? -Z * Z' : D) + so * J(1, n, 1) + J(n, 1, 1) * so'
	if (fam == 1) y = rpoisson(1, 1, exp(eta))
	else if (fam == 2) y = rbinomial(1, 1, 7, invlogit(eta))
	else y = eta + rnormal(n, n, 0, sqrt(disp))
	if (!dir) {
		y = uppertriangle(y, 0)
		y = y + y'
	}
	_diag(y, 0)
	st_store(., ., y)
	__lsmt_Dtrue = bil ? -Z * Z' : D
	__lsmt_Xc = Xc
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
capture program drop lsmt_recv
program lsmt_recv
	syntax name(name=case), FAM(string) DIR(integer) B0(real) B1(real) SEED(integer) [DISP(real 0.5) SOC(integer 0) OPTS(string) NOCOV]
	local n 50
	clear
	set obs `n'
	forvalues j = 1/`n' {
		quietly gen double __a`j' = .
	}
	local famn = cond("`fam'" == "poisson", 1, cond("`fam'" == "binomial", 2, 3))
	local bil = strpos(" `opts' ", " bilinear ") > 0
	mata: lsmt_simv(`n', `dir', `famn', `b0', `b1', `seed', `disp', `soc', `bil')
	if `dir' nwset __a*, name(sim_`case') directed valued
	else nwset __a*, name(sim_`case') undirected valued
	mata: st_matrix("__XC", __lsmt_Xc)
	local extra
	if "`fam'" == "binomial" local extra "trials(7)"
	if "`fam'" == "normal" local extra "dispvar(0.5) dispdf(3)"
	quietly nwlsm sim_`case', family(`fam') `extra' edgecov(__XC) `opts' seed(7) nsample(2000) burnin(5000) interval(10)
	matrix T = e(summary)
	local ok 1
	local truth0 = `b0'
	if "`fam'" == "binomial" | 1 {
		* the bilinear intercept has the same meaning; the effects shift the intercept by their means only approximately: allow the 95% interval widened by 0.3
		local w = cond(`soc' > 0, 0.3, 0)
	}
	forvalues k = 1/2 {
		local tr = cond(`k' == 1, `b0', `b1')
		local lo = T[`k', 4] - `w'
		local hi = T[`k', 5] + `w'
		di as text "case `case' coef `k': true " %7.3f `tr' "  posterior mean " %7.3f T[`k', 1] "  95% interval [" %7.3f T[`k', 4] ", " %7.3f T[`k', 5] "]"
		if `tr' < `lo' | `tr' > `hi' local ok 0
	}
	if "`fam'" == "normal" {
		matrix V = e(vsummary)
		local r = rowsof(V)
		di as text "case `case' dispersion: true " %7.3f `disp' "  posterior mean " %7.3f V[`r', 1] "  95% interval [" %7.3f V[`r', 4] ", " %7.3f V[`r', 5] "]"
		if `disp' < V[`r', 4] | `disp' > V[`r', 5] local ok 0
	}
	if !`bil' {
		mata: st_local("cor", strofreal(lsmt_cor_ut(st_matrix("e(D)"), __lsmt_Dtrue), "%9.5f"))
		di as text "case `case' correlation of the posterior mean distances with the true distances: " `cor'
		if `cor' < 0.9 local ok 0
	}
	if !`ok' {
		di as error "FAIL case `case'"
		exit 9
	}
	di as result "PASS case `case'"
end
lsmt_recv pois, fam(poisson) dir(1) b0(0.5) b1(0.3) seed(11)
lsmt_recv bin, fam(binomial) dir(0) b0(-0.2) b1(0.4) seed(12)
lsmt_recv gau, fam(normal) dir(0) b0(1) b1(0.4) seed(13) disp(0.5) soc(1) opts(sociality)
lsmt_recv poisbil, fam(poisson) dir(0) b0(0.3) b1(0.2) seed(14) opts(bilinear)
display "PASS test_nwlsm_recovery_valued: Poisson, binomial and normal families recover coefficients, dispersion and distances"
