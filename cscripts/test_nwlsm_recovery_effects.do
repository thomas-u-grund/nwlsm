* Recovery with random effects: networks simulated from known coefficients, positions and sender / receiver / sociality effects (true standard
* deviations 0.8 and 0.6). The 95% interval of the covariate coefficient must contain the true value, the posterior mean effects must correlate with
* the true effects (at least 0.7: an effect rests on the ties of one actor only, and the true standard deviations are small), the 95% intervals of the variances must contain the true variances, the distance matrix must be recovered
* (correlation at least 0.75 for the distances, 0.8 for the inner products), and the bilinear model must recover the inner products.
do cscripts/nwlsm_testlib.do
mata:
mata set matastrict off
void lsmt_simfx(real scalar n, real scalar dir, real scalar d, real scalar seed, real scalar kind)
{
	real matrix Z, D, Xc, eta, y
	real colvector x, se, re
	real scalar i
	external real matrix __lsmt_Dtrue, __lsmt_Xc, __lsmt_Itrue
	external real colvector __lsmt_se, __lsmt_re

	rseed(seed)
	Z = rnormal(n, d, 0, 1.6)
	x = rnormal(n, 1, 0, 1)
	se = rnormal(n, 1, 0, 0.8)
	re = rnormal(n, 1, 0, 0.6)
	D = J(n, n, 0)
	for (i = 1; i <= d; i++) D = D + (Z[., i] * J(1, n, 1) - J(n, 1, 1) * Z[., i]') :^ 2
	D = sqrt(D)
	Xc = abs(x * J(1, n, 1) - J(n, 1, 1) * x')
	if (kind == 1) eta = 0.5 :+ 0.7 * Xc - D + se * J(1, n, 1) + J(n, 1, 1) * re'
	else if (kind == 2) eta = 0.5 :+ 0.7 * Xc - D + se * J(1, n, 1) + J(n, 1, 1) * se'
	else eta = 0.5 :+ 0.7 * Xc + 0.4 * (Z * Z')
	y = (runiform(n, n) :< invlogit(eta))
	if (!dir) {
		y = uppertriangle(y, 0)
		y = y + y'
	}
	_diag(y, 0)
	st_store(., ., y)
	__lsmt_Dtrue = D
	__lsmt_Itrue = Z * Z'
	__lsmt_Xc = Xc
	__lsmt_se = se
	__lsmt_re = re
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
capture program drop lsmt_recfx
program lsmt_recfx
	syntax, N(integer) DIR(integer) D(integer) SEED(integer) KIND(integer) OPT(string) [BIL]
	clear
	set obs `n'
	forvalues j = 1/`n' {
		quietly gen double __a`j' = .
	}
	mata: lsmt_simfx(`n', `dir', `d', `seed', `kind')
	if `dir' nwset __a*, name(rec) directed replace
	else nwset __a*, name(rec) undirected replace
	mata: st_matrix("__XC", __lsmt_Xc)
	quietly nwlsm rec, dim(`d') edgecov(__XC) `opt' seed(7) nsample(3000) burnin(8000)
	local ok 1
	matrix T = e(summary)
	di as text "  kind `kind' n=`n': covariate coefficient: true 0.7, posterior mean " %7.4f T[2, 1] "  95% interval [" %7.4f T[2, 4] ", " %7.4f T[2, 5] "]"
	if 0.7 < T[2, 4] | 0.7 > T[2, 5] local ok 0
	matrix V = e(vsummary)
	local nv = rowsof(V)
	local rn : rownames V
	forvalues k = 2/`nv' {
		local nmv : word `k' of `rn'
		local tv = cond("`nmv'" == "sender_var", 0.64, cond("`nmv'" == "receiver_var", 0.36, 0.64))
		di as text "  `nmv': true " %5.2f `tv' "  posterior mean " %7.4f V[`k', 1] "  95% interval [" %7.4f V[`k', 4] ", " %7.4f V[`k', 5] "]"
		if `tv' < V[`k', 4] | `tv' > V[`k', 5] local ok 0
	}
	if "`bil'" == "" {
		mata: st_local("cor", strofreal(lsmt_cor_ut(st_matrix("e(D)"), __lsmt_Dtrue), "%8.4f"))
		di as text "  correlation of true and posterior mean distances: " `cor'
		if `cor' < 0.75 local ok 0
		if strpos("`opt'", "sender") {
			matrix S = e(sender)
			mata: st_local("cs", strofreal(correlation((st_matrix("S")[., 1], __lsmt_se))[1, 2], "%8.4f"))
			matrix R = e(receiver)
			mata: st_local("cr", strofreal(correlation((st_matrix("R")[., 1], __lsmt_re))[1, 2], "%8.4f"))
			di as text "  correlation of true and posterior mean sender effects: " `cs' ", receiver effects: " `cr'
			if `cs' < 0.7 | `cr' < 0.7 local ok 0
		}
		else {
			matrix S = e(sociality)
			mata: st_local("cs", strofreal(correlation((st_matrix("S")[., 1], __lsmt_se))[1, 2], "%8.4f"))
			di as text "  correlation of true and posterior mean sociality effects: " `cs'
			if `cs' < 0.7 local ok 0
		}
	}
	else {
		mata: st_local("cor", strofreal(lsmt_cor_ut(__lsmt_ip(st_matrix("e(Z)")), __lsmt_Itrue), "%8.4f"))
		di as text "  correlation of true and estimated inner products (MKL positions): " `cor'
		if `cor' < 0.8 local ok 0
	}
	if !`ok' {
		di as error "FAIL recovery with effects kind `kind'"
		exit 9
	}
	di as result "PASS recovery with effects kind `kind' n=`n'"
end
mata:
real matrix __lsmt_ip(real matrix Z) return(Z * Z')
end
lsmt_recfx, n(70) dir(1) d(2) seed(11) kind(1) opt(sender receiver)
lsmt_recfx, n(70) dir(0) d(2) seed(12) kind(2) opt(sociality)
lsmt_recfx, n(70) dir(0) d(2) seed(13) kind(3) opt(bilinear) bil
display "PASS test_nwlsm_recovery_effects"
