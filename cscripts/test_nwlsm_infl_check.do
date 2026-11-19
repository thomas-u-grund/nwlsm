* Influence model, update check: the plugin's outcome log likelihood, the weighted neighbour means m_k and the incremental change of the outcome
* log likelihood when one actor moves (the quantity added to the Metropolis acceptance ratio, including the shortcuts for rows whose neighbour set
* does not change) equal a direct computation in Mata that selects the G nearest actors by sorting. Random states: n = 5 to 40, inverse distance and
* exponential weights, all or G nearest actors, 1 to 4 extra regressors, small and large moves, directed and undirected ties.
do cscripts/nwlsm_testlib.do
set seed 20261120
mata:
mata set matastrict off
// outcome log likelihood and the means m_k by direct computation
real matrix lsmt_infl(real matrix Z, real colvector yo, real colvector yl, real matrix Xo, real colvector ob, real scalar s2, real scalar G, real scalar wt)
{
	real scalar n, k, j, Gm, num, den, f
	real matrix D
	real colvector m, o, res, w
	n = rows(Z)
	D = lsmt_pd(Z)
	Gm = (G <= 0 | G >= n - 1) ? n - 1 : G
	m = J(n, 1, 0)
	for (k = 1; k <= n; k++) {
		o = order(D[k, .]', 1)
		o = select(o, o :!= k)
		o = o[1..Gm]
		w = wt == 0 ? 1 :/ D[k, o]' : exp(-D[k, o]')
		m[k] = sum(w :* yl[o]) / sum(w)
	}
	res = yo - Xo * ob[1..(cols(Xo))] - ob[cols(Xo) + 1] * m
	return((-0.5 * n * ln(2 * pi() * s2) - sum(res :^ 2) / (2 * s2) \ m))
}
void lsmt_run()
{
	maxe = 0; ntest = 0
	for (rep = 1; rep <= 120; rep++) {
		n = 5 + floor(runiform(1, 1) * 36)
		dir = mod(rep, 2)
		d = 1 + mod(rep, 3)
		y = (runiform(n, n) :< 0.3)
		if (!dir) {
			y = uppertriangle(y, 0); y = y + y'
		}
		_diag(y, 0)
		p = 2 + mod(rep, 4)
		Xo = (J(n, 1, 1), rnormal(n, p - 1, 0, 1))
		yl = Xo[., 2]
		yo = rnormal(n, 1, 0, 2)
		ob = rnormal(p + 1, 1, 0, 1)
		s2 = 0.3 + runiform(1, 1) * 2
		wt = mod(floor(rep / 2), 2)
		G = (mod(rep, 3) == 0) ? 0 : 1 + floor(runiform(1, 1) * (n - 1))
		Z = rnormal(n, d, 0, 1.5)
		i = 1 + floor(runiform(1, 1) * n)
		znew = (mod(rep, 4) < 2) ? Z[i, .] + rnormal(1, d, 0, 0.1) : rnormal(1, d, 0, 1.5)
		io = (yo, yl, Xo)
		isv = (G, wt, 100, 0.01, 0.01, s2, ob')
		r = LsmCheckInfl(y, dir, J(0, n, .), d, 0, Z, 1, io, isv, i, znew)
		a0 = lsmt_infl(Z, yo, yl, Xo, ob, s2, G, wt)
		Z2 = Z
		Z2[i, .] = znew
		a1 = lsmt_infl(Z2, yo, yl, Xo, ob, s2, G, wt)
		dd = a1[1] - a0[1]
		e = (abs(r[1] - a0[1]), max(abs(r[2..(n + 1)]' - a0[2..(n + 1)])), abs(r[n + 2] - dd), abs(r[n + 3] - dd), abs(r[n + 4] - dd))
		maxe = max((maxe, e))
		ntest++
		if (max(e) > 1e-7 * (1 + abs(a0[1]))) {
			printf("MISMATCH rep %f n %f G %f wt %f dir %f: errors %s\n", rep, n, G, wt, dir, invtokens(strofreal(e, "%9.3g")))
			exit(9)
		}
	}
	printf("PASS: %f random states, largest error %g (outcome log likelihood, weighted means, incremental change, recomputed change, state after commit)\n", ntest, maxe)
}
end
mata: lsmt_run()
