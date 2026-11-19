* Likelihood and single-update check: the plugin's log likelihood equals a direct recomputation in Mata, and the incremental
* log-likelihood change of a one-actor update (the quantity in the Metropolis acceptance ratio) equals the difference of two
* full recomputations. Random states, directed and undirected, d = 1, 2, 3, covariates, missing ties.
do cscripts/nwlsm_testlib.do
set seed 20261119
mata:
mata set matastrict off
void lsmt_run()
{
	t_maxerr1 = 0; t_maxerr2 = 0; ntests = 0
	for (rep = 1; rep <= 60; rep++) {
		n = 4 + floor(runiform(1, 1) * 9)
		dir = mod(rep, 2)
		d = 1 + mod(rep, 3)
		nc = mod(floor(rep / 3), 3)
		y = (runiform(n, n) :< 0.35)
		mm = (runiform(n, n) :< 0.1) :* (rep > 40)
		if (!dir) {
			y = uppertriangle(y, 0); y = y + y'
			mm = uppertriangle(mm, 0); mm = mm + mm'
		}
		y = y :/ (!mm)
		_diag(y, 0)
		X = J(0, n, .)
		for (k = 1; k <= nc; k++) {
			Xk = rnormal(n, n, 0, 1)
			if (!dir) {
				Xk = uppertriangle(Xk, 0)
				Xk = Xk + Xk'
			}
			_diag(Xk, 0)
			X = X \ Xk
		}
		beta = rnormal(1, nc + 1, 0, 1)
		Z = rnormal(n, d, 0, 1.2)
		i = 1 + floor(runiform(1, 1) * n)
		znew = rnormal(1, d, 0, 1.5)
		s2 = 0.5 + runiform(1, 1) * 3
		r = LsmCheck(y, dir, X, d, beta, Z, s2, i, znew)
		ll0 = lsmt_ll(y, dir, X, beta, Z)
		Z2 = Z
		Z2[i, .] = znew
		dll = lsmt_ll(y, dir, X, beta, Z2) - ll0
		lp = -sum(Z :^ 2) / (2 * s2) - 0.5 * n * d * ln(2 * pi() * s2)
		e1 = abs(r[1] - ll0); e2 = abs(r[3] - dll); e3 = abs(r[2] - lp)
		if (max((e1, e2, e3)) > 1e-8) {
			printf("MISMATCH rep %f n %f dir %f d %f nc %f: ll %g %g, dll %g %g, prior %g %g\n", rep, n, dir, d, nc, r[1], ll0, r[3], dll, r[2], lp)
			exit(9)
		}
		t_maxerr1 = max((t_maxerr1, e1)); t_maxerr2 = max((t_maxerr2, e2, e3))
		ntests++
	}
	printf("cases: %f, max |ll error| %g, max |dll or log prior error| %g\n", ntests, t_maxerr1, t_maxerr2)
}
lsmt_run()
end
display "PASS test_nwlsm_check: log likelihood, one-actor update and prior density agree with the direct computation (60 random states)"
