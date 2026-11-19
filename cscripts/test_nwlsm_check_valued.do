* Valued families: the plugin's log likelihood (Poisson, binomial, normal), the incremental changes of one-actor and one-effect updates and the analytic
* gradient of the minimum-KL objective equal direct recomputations / central differences. 90 random states (directed and undirected, d = 1 to 3,
* covariates, missing values, effects, Euclidean and bilinear).
do cscripts/nwlsm_testlib.do
set seed 20261120
mata:
mata set matastrict off
real scalar lsmt_llv(real matrix y, real scalar dir, real matrix X, real rowvector beta, real matrix Z, real matrix E, real rowvector opt, real scalar fam, real matrix tr, real scalar disp)
{
	real scalar n, k, i, j, eta, t, lt

	n = rows(y)
	t = 0
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j | (!dir & j < i)) continue
		if (y[i, j] >= .) continue
		eta = beta[1]
		for (k = 1; k <= rows(X) / n; k++) eta = eta + beta[k + 1] * X[(k - 1) * n + i, j]
		if (opt[1]) eta = eta + E[i, 1]
		if (opt[2]) eta = eta + E[j, 2]
		if (opt[3]) eta = eta + E[i, 3] + E[j, 3]
		if (opt[4]) lt = -sum(Z[i, .] :* Z[j, .])
		else lt = sqrt(sum((Z[i, .] - Z[j, .]) :^ 2))
		eta = eta - lt
		if (fam == 1) t = t + y[i, j] * eta - exp(eta) - lngamma(y[i, j] + 1)
		else if (fam == 2) t = t + y[i, j] * eta - tr[i, j] * ln(1 + exp(eta)) + lngamma(tr[i, j] + 1) - lngamma(y[i, j] + 1) - lngamma(tr[i, j] - y[i, j] + 1)
		else t = t - (y[i, j] - eta)^2 / (2 * disp) - 0.5 * ln(2 * pi() * disp)
	}
	return(t)
}
void lsmt_runv()
{
	maxe = 0; maxg = 0; ntests = 0
	for (rep = 1; rep <= 90; rep++) {
		fam = 1 + mod(rep, 3)
		n = 4 + floor(runiform(1, 1) * 9)
		dir = mod(floor(rep / 3), 2)
		d = 1 + mod(rep, 3)
		nc = mod(floor(rep / 3), 3)
		tr = 2 :+ floor(runiform(n, n) * 6)
		if (!dir) {
			tr = uppertriangle(tr, 0)
			tr = tr + tr'
			tr = tr :+ (tr :== 0)
		}
		if (fam == 1) y = floor(runiform(n, n) * 4)
		else if (fam == 2) y = floor(runiform(n, n) :* (tr :+ 1))
		else y = rnormal(n, n, 0.5, 1.3)
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
		beta = rnormal(1, nc + 1, 0, 0.5)
		Z = rnormal(n, d, 0, 1.0)
		i = 1 + floor(runiform(1, 1) * n)
		znew = rnormal(1, d, 0, 1.0)
		s2 = 0.5 + runiform(1, 1) * 3
		disp = 0.4 + runiform(1, 1) * 2
		opt = (0, 0, 0, 0, 1, 3, 1, 3, 1, 3)
		if (rep > 30) {
			if (dir) opt[1..3] = (runiform(1, 3) :< 0.6)
			else opt[3] = (runiform(1, 1) < 0.7)
		}
		opt[4] = (mod(rep, 5) == 0)
		E = rnormal(n, 3, 0, 0.4) :* (opt[1..3] :!= 0)
		t = 0
		if (sum(opt[1..3]) > 0) {
			t = selectindex(opt[1..3])
			t = t[1 + floor(runiform(1, 1) * cols(t))]
		}
		delta = rnormal(1, 1, 0, 0.5)
		fm = (fam, 1, 3, (fam == 3 ? disp : 0))
		r = LsmCheck(y, dir, X, d, beta, Z, s2, i, znew, opt, E, t, delta, fm, tr)
		ll0 = lsmt_llv(y, dir, X, beta, Z, E, opt, fam, tr, disp)
		Z2 = Z
		Z2[i, .] = znew
		dll = lsmt_llv(y, dir, X, beta, Z2, E, opt, fam, tr, disp) - ll0
		e1 = abs(r[1] - ll0) / (1 + abs(ll0)); e2 = abs(r[3] - dll) / (1 + abs(dll)); e4 = 0
		if (t > 0) {
			E2 = E
			E2[i, t] = E2[i, t] + delta
			dle = lsmt_llv(y, dir, X, beta, Z, E2, opt, fam, tr, disp) - ll0
			e4 = abs(r[4] - dle) / (1 + abs(dle))
		}
		if (max((e1, e2, e4)) > 1e-9) {
			printf("MISMATCH rep %f fam %f n %f dir %f d %f nc %f: ll %g %g, dll %g %g, effect err %g\n", rep, fam, n, dir, d, nc, r[1], ll0, r[3], dll, e4)
			exit(9)
		}
		g = LsmCheckMKL(y, dir, X, d, beta, Z, opt, E, fm, tr)
		P = (cols(g) - 1) / 2
		eg = max(abs(g[2..(1 + P)] - g[(2 + P)..(1 + 2 * P)]) :/ (1 :+ abs(g[(2 + P)..(1 + 2 * P)])))
		if (eg > 1e-6) {
			printf("GRADIENT MISMATCH rep %f fam %f: %g\n", rep, fam, eg)
			exit(9)
		}
		maxg = max((maxg, eg)); maxe = max((maxe, e1, e2, e4)); ntests++
	}
	printf("cases: %f, max relative error of ll and updates %g, max relative gradient error %g\n", ntests, maxe, maxg)
}
lsmt_runv()
end
display "PASS test_nwlsm_check_valued: Poisson, binomial and normal likelihoods, incremental updates and MKL gradients agree with the direct computation (90 random states)"
