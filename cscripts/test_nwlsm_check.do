* Likelihood and single-update check: the plugin's log likelihood equals a direct recomputation in Mata, and the incremental
* log-likelihood changes of a one-actor update and of a one-effect update (the quantities in the Metropolis acceptance ratios)
* equal the differences of two full recomputations; the analytic gradient of the MKL objective equals a central difference.
* Random states, directed and undirected, d = 1, 2, 3, covariates, missing ties, sender/receiver/sociality effects, euclidean and bilinear.
do cscripts/nwlsm_testlib.do
set seed 20261119
mata:
mata set matastrict off
void lsmt_run()
{
	t_maxerr1 = 0; t_maxerr2 = 0; t_maxerr3 = 0; t_maxg = 0; ntests = 0; neff = 0; nbil = 0
	for (rep = 1; rep <= 90; rep++) {
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
		// effects: in a directed network any mix of sender, receiver, sociality; undirected: sociality only
		opt = (0, 0, 0, 0, 1, 3, 1, 3, 1, 3)
		if (rep > 30) {
			if (dir) opt[1..3] = (runiform(1, 3) :< 0.6)
			else opt[3] = (runiform(1, 1) < 0.7)
		}
		opt[4] = (mod(rep, 5) == 0)
		E = rnormal(n, 3, 0, 0.7) :* (opt[1..3] :!= 0)
		t = 0
		if (sum(opt[1..3]) > 0) {
			t = selectindex(opt[1..3])
			t = t[1 + floor(runiform(1, 1) * cols(t))]
		}
		delta = rnormal(1, 1, 0, 0.8)
		r = LsmCheck(y, dir, X, d, beta, Z, s2, i, znew, opt, E, t, delta)
		ll0 = lsmt_ll(y, dir, X, beta, Z, E, opt)
		Z2 = Z
		Z2[i, .] = znew
		dll = lsmt_ll(y, dir, X, beta, Z2, E, opt) - ll0
		lp = -sum(Z :^ 2) / (2 * s2) - 0.5 * n * d * ln(2 * pi() * s2)
		e1 = abs(r[1] - ll0); e2 = abs(r[3] - dll); e3 = abs(r[2] - lp)
		e4 = 0
		if (t > 0) {
			E2 = E
			E2[i, t] = E2[i, t] + delta
			e4 = abs(r[4] - (lsmt_ll(y, dir, X, beta, Z, E2, opt) - ll0))
			neff++
		}
		if (max((e1, e2, e3, e4)) > 1e-8) {
			printf("MISMATCH rep %f n %f dir %f d %f nc %f opt %s: ll %g %g, dll %g %g, prior %g %g, effect dll err %g\n", rep, n, dir, d, nc, invtokens(strofreal(opt[1..4])), r[1], ll0, r[3], dll, r[2], lp, e4)
			exit(9)
		}
		nbil = nbil + opt[4]
		// gradient of the MKL objective
		g = LsmCheckMKL(y, dir, X, d, beta, Z, opt, E)
		P = (cols(g) - 1) / 2
		ga = g[2..(1 + P)]; gn = g[(2 + P)..(1 + 2 * P)]
		eg = max(abs(ga - gn) :/ (1 :+ abs(gn)))
		if (eg > 1e-6) {
			printf("GRADIENT MISMATCH rep %f: %g\n", rep, eg)
			exit(9)
		}
		t_maxg = max((t_maxg, eg))
		t_maxerr1 = max((t_maxerr1, e1)); t_maxerr2 = max((t_maxerr2, e2, e3, e4))
		ntests++
	}
	printf("cases: %f (effect updates %f, bilinear %f), max |ll error| %g, max |dll, log prior or effect dll error| %g, max relative gradient error %g\n", ntests, neff, nbil, t_maxerr1, t_maxerr2, t_maxg)
}
lsmt_run()
end
display "PASS test_nwlsm_check: log likelihood, one-actor and one-effect updates, prior density and MKL gradient agree with the direct computation (90 random states)"
