* Recovery: networks simulated from known coefficients and positions; the fit must cover the true coefficients with the 95% interval and
* reproduce the true distance matrix (correlation of the upper triangles). Undirected n = 60 with a node covariate (absdiff and edgecov
* must give identical results), directed n = 50 with one edge covariate, d = 2 each; and a d = 3 case.
do cscripts/nwlsm_testlib.do
capture program drop lsmt_sim
* lsmt_sim: simulates y (stored as variables __a1..__an of a fresh dataset) and the true distance matrix, covariate matrix and x in Mata externals
program lsmt_sim
	syntax, N(integer) DIR(integer) D(integer) B0(real) B1(real) SEED(integer) ZSD(real)
	clear
	set obs `n'
	forvalues j = 1/`n' {
		quietly gen double __a`j' = .
	}
	mata: lsmt_simulate(`n', `dir', `d', `b0', `b1', `seed', `zsd')
end
mata:
void lsmt_simulate(real scalar n, real scalar dir, real scalar d, real scalar b0, real scalar b1, real scalar seed, real scalar zsd)
{
	real matrix Z, D, Xc, eta, y
	real colvector x
	real scalar i, j
	external real matrix __lsmt_Dtrue, __lsmt_Xc
	external real colvector __lsmt_x

	rseed(seed)
	Z = rnormal(n, d, 0, zsd)
	x = rnormal(n, 1, 0, 1)
	D = J(n, n, 0)
	for (i = 1; i <= d; i++) D = D + (Z[., i] * J(1, n, 1) - J(n, 1, 1) * Z[., i]') :^ 2
	D = sqrt(D)
	Xc = abs(x * J(1, n, 1) - J(n, 1, 1) * x')
	if (!dir) Xc = Xc
	eta = b0 :+ b1 * Xc - D
	y = (runiform(n, n) :< invlogit(eta))
	if (!dir) {
		y = uppertriangle(y, 0)
		y = y + y'
	}
	_diag(y, 0)
	st_store(., ., y)
	__lsmt_Dtrue = D
	__lsmt_Xc = Xc
	__lsmt_x = x
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
capture program drop lsmt_rec
program lsmt_rec
	syntax, N(integer) DIR(integer) D(integer) B0(real) B1(real) SEED(integer) ZSD(real) [ABSdiff]
	lsmt_sim, n(`n') dir(`dir') d(`d') b0(`b0') b1(`b1') seed(`seed') zsd(`zsd')
	if `dir' nwset __a*, name(rec) directed replace
	else nwset __a*, name(rec) undirected replace
	mata: st_matrix("__XC", __lsmt_Xc)
	quietly nwlsm rec, dim(`d') edgecov(__XC) seed(7) nsample(2000) burnin(5000)
	local ok 1
	matrix b = e(b)
	matrix T = e(summary)
	local k = 0
	foreach tb in `b0' `b1' {
		local ++k
		local lo = T[`k', 4]
		local hi = T[`k', 5]
		di as text "  n=`n' dir=`dir' d=`d' coef `k': true " %6.3f `tb' "  posterior mean " %7.4f b[1, `k'] "  95% interval [" %7.4f `lo' ", " %7.4f `hi' "]"
		if `tb' < `lo' | `tb' > `hi' local ok 0
	}
	mata: st_local("cor", strofreal(lsmt_cor_ut(st_matrix("e(D)"), __lsmt_Dtrue), "%8.4f"))
	di as text "  correlation of true and posterior mean distances: " `cor'
	if `cor' < cond(`d' >= 3, 0.8, 0.9) local ok 0
	if !`ok' {
		di as error "FAIL recovery n=`n' dir=`dir' d=`d'"
		exit 9
	}
	if "`absdiff'" != "" {
		* node-attribute form must reproduce the edgecov fit exactly
		matrix bb = e(b)
		preserve
		clear
		set obs `n'
		mata: st_addvar("double", "xx"); st_store(., "xx", __lsmt_x)
		quietly nwlsm rec, dim(`d') absdiff(xx) seed(7) nsample(2000) burnin(5000)
		matrix b2 = e(b)
		restore
		assert reldif(bb[1, 1], b2[1, 1]) < 1e-12 & reldif(bb[1, 2], b2[1, 2]) < 1e-12
		di as result "  absdiff(xx) gives the same fit as the equivalent edgecov()"
	}
	di as result "PASS recovery n=`n' dir=`dir' d=`d'"
end
lsmt_rec, n(60) dir(0) d(2) b0(2) b1(0.8) seed(101) zsd(2) absdiff
lsmt_rec, n(50) dir(1) d(2) b0(2) b1(-0.6) seed(202) zsd(1.8)
lsmt_rec, n(60) dir(0) d(3) b0(2) b1(0.5) seed(303) zsd(2)
display "PASS test_nwlsm_recovery"
