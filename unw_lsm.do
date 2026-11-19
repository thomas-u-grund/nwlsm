/*
 * Copyright (C) 2026 Thomas U. Grund
 * This file is part of nwlsm, a Stata package for latent space models for networks. It follows the latent space
 * model of Hoff, Raftery and Handcock (2002) as implemented in the R package latentnet (Krivitsky and Handcock) and
 * is written independently.
 * nwlsm is free software: you can redistribute it and/or modify it under the terms of the GNU General
 * Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
 */
/*
	unw_lsm.do -- Mata glue of nwlsm (built into lib/lnwlsm.mlib). The sampler runs in the compiled plugin native/lsm_sim.c;
	Mata prepares the data and the start values (classical multidimensional scaling of the geodesic distances, a logit for the
	coefficients), writes the model into the plugin's frame, and summarises the posterior sample (coefficients, mean distances,
	Procrustes-aligned positions). There is no Mata fallback: without a plugin the command stops (r(601), r(3495)).
*/
capture mata: mata drop LsmPlatform()
capture mata: mata drop LsmPluginFilename()
capture mata: mata drop LsmPluginPath()
capture mata: mata drop LsmNativeAvailable()
capture mata: mata drop LSM_PROTOCOL_NEEDED()
capture mata: mata drop LsmNat()
capture mata: mata drop lsm_geodesic()
capture mata: mata drop lsm_mds()
capture mata: mata drop lsm_start()
capture mata: mata drop LsmGof()
capture mata: mata drop LsmTraceMat()
capture mata: mata drop LsmPPCheck()
capture mata: mata drop lsm_gofsum()
capture mata: mata drop lsm_vstart()
capture mata: mata drop lsm_checkvalues()
capture mata: mata drop lsm_priorvar()
capture mata: mata drop lsm_mcse()
capture mata: mata drop lsm_rhat()
capture mata: mata drop lsm_perms()
capture mata: mata drop lsm_relabel()
capture mata: mata drop lsm_procR()
capture mata: mata drop lsm_procrustes()
capture mata: mata drop lsm_ymat()
capture mata: mata drop lsm_node_covs()
capture mata: mata drop LsmFit()
capture mata: mata drop LsmCheck()
capture mata: mata drop LsmSample()
capture mata: mata drop LsmSaveDraws()
capture mata: mata drop lsm_ess()
capture mata: mata drop lsm_geweke()
capture mata: mata drop LsmDiag()
capture mata: mata drop lsm_vartable()
capture mata: mata drop LsmCheckMKL()
capture mata: mata drop LsmCheckInfl()
capture mata: mata drop lsm_ztrans()
capture mata: mata drop LsmPlotPrep()
capture mata: mata drop LsmPlotStore()
capture mata: mata drop LsmBic()
capture mata: mata drop LsmBicFit()
capture mata: mata drop lsm_condmle()
capture mata: mata drop lsm_kmeans()
capture mata: mata drop lsm_em_vii()

mata:
mata set matastrict off

string scalar LsmPlatform()
{
	string scalar os

	os = st_global("c(os)")
	if (os == "Windows") return("windows")
	if (os == "MacOSX") return("macos")
	if (substr(st_global("c(machine_type)"), 1, 3) == "Mac") return("macos")
	return("unix")
}

string scalar LsmPluginFilename() return("lsm_sim_" + LsmPlatform() + ".plugin")

// (1) findfile() on the file name alone (after a net install), (2) lib/plugins/<os>/ next to nwlsm.ado (a checkout)
string scalar LsmPluginPath()
{
	string scalar fname, found, full, dir, fn

	fname = LsmPluginFilename()
	found = findfile(fname)
	if (found != "") return(found)
	full = findfile("nwlsm.ado")
	if (full == "") return("")
	pathsplit(full, dir, fn)
	return(pathjoin(pathjoin(dir, "lib"), pathjoin("plugins", pathjoin(LsmPlatform(), fname))))
}

real scalar LsmNativeAvailable()
{
	string scalar p

	p = LsmPluginPath()
	if (p == "") return(0)
	return(fileexists(p))
}

real scalar LSM_PROTOCOL_NEEDED() return(6)

class LsmNat {
	real scalar isopen, L
	string scalar origframe
	void open()
	void close()
	void setup()
	real colvector call()
	void append()
}

void LsmNat::open()
{
	string scalar p, fname, q
	real scalar rc, ver

	isopen = 0
	fname = LsmPluginFilename()
	p = LsmPluginPath()
	if (p == "" | !fileexists(p)) {
		errprintf("{p 0 0 2}nwlsm: the compiled plugin {bf:%s} was not found. nwlsm fits every model with this plugin and has no Mata fallback; reinstall the package with {bf:net install nwlsm, replace}.{p_end}\n", fname)
		exit(601)
	}
	origframe = st_framecurrent()
	stata("capture frame drop __nwlsm_native")
	stata("frame create __nwlsm_native")
	st_framecurrent("__nwlsm_native")
	isopen = 1
	st_addobs(4)
	(void) st_addvar("double", "v1")
	(void) st_addvar("double", "v2")
	q = char(34)
	stata("capture program drop lsmnative")
	stata("capture program lsmnative, plugin using(" + q + p + q + ")")
	rc = _stata("plugin call lsmnative v1 v2, " + q + "99 0" + q, 1)
	ver = .
	if (rc == 0) ver = st_data(1, 2)
	if (rc != 0 | ver != LSM_PROTOCOL_NEEDED()) {
		close()
		errprintf("{p 0 0 2}nwlsm: the installed plugin {bf:%s} does not match this version of nwlsm (or it cannot be loaded). Reinstall the package with {bf:net install nwlsm, replace}.{p_end}\n", fname)
		exit(3495)
	}
}

void LsmNat::close()
{
	if (isopen) {
		st_framecurrent(origframe)
		stata("capture frame drop __nwlsm_native")
		isopen = 0
	}
}

// writes more numbers behind the model block (read by the goodness-of-fit call)
void LsmNat::append(real colvector v)
{
	if (L + rows(v) > st_nobs()) st_addobs(L + rows(v) - st_nobs())
	st_store((L + 1::L + rows(v)), 1, v)
}

// y: n x n (0/1, missing = .), X: (nc*n) x n, the covariate matrices stacked, d: dimension, prior: (scale, df, s2 start),
// pvar, b0: 1 x K, Z0: n x d, flags: (sender, receiver, sociality, bilinear), eprior: 3 x 3 (scale, df, start variance of the effects),
// E0: n x 3 start values of the effects
void LsmNat::setup(real matrix y, real scalar directed, real matrix X, real scalar d, real rowvector prior,
		real rowvector pvar, real rowvector b0, real matrix Z0, real rowvector flags, real matrix eprior, real matrix E0, real rowvector cl, | real rowvector fm, real matrix tr, real matrix io, real rowvector isv)
{
	real scalar n, nc, k, ion, p
	real colvector blob, yb

	n = rows(y)
	nc = rows(X) / n
	if (args() < 13) fm = (0, 0, 0, 0)
	if (args() < 14) tr = J(0, 0, .)
	if (args() < 16) {
		io = J(0, 0, .)
		isv = J(1, 6, 0)
	}
	ion = (rows(io) > 0)
	p = ion ? cols(io) - 2 : 0
	// y block: binary: 0 / 1 / -1 (missing); valued families: the values (0 for missing), the observed indicators, the trials (binomial)
	if (fm[1] == 0) yb = vec(editmissing(y, -1)')
	else {
		yb = vec(editmissing(y, 0)') \ vec((y :< .)')
		if (fm[1] == 2) yb = yb \ vec(tr')
	}
	blob = (LSM_PROTOCOL_NEEDED() \ n \ directed \ nc \ d \ (flags[1] + 2 * flags[2] + 4 * flags[3]) \ flags[4] \ cl[1] \ cl[2] \ cl[3] \ cl[4] \ cl[5] \ fm[1] \ fm[2] \ fm[3] \ fm[4] \ ion \ (ion ? isv[1] : 0) \ (ion ? isv[2] : 0) \ p) \ prior' \ pvar' \ b0' \ vec(Z0') \ vec(eprior') \ vec(E0) \ yb
	if (nc > 0) {
		for (k = 1; k <= nc; k++) blob = blob \ vec(editmissing(X[((k - 1) * n + 1)..(k * n), .], 0)')
	}
	// influence model: prior variance of the coefficients, inverse gamma (a, b) of the error variance, start of the error variance, start of the coefficients (p + 1),
	// outcome, lagged outcome, the p regressors (by column)
	if (ion) blob = blob \ isv[3]' \ isv[4]' \ isv[5]' \ isv[6]' \ isv[7..(7 + p)]' \ io[., 1] \ io[., 2] \ vec(io[., 3..(2 + p)])
	L = rows(blob)
	if (L > st_nobs()) st_addobs(L - st_nobs())
	st_store((1::L), 1, blob)
}

real colvector LsmNat::call(string scalar arg, real scalar nrows)
{
	real scalar need

	need = max((L, nrows, 4))
	if (need > st_nobs()) st_addobs(need - st_nobs())
	stata("plugin call lsmnative v1 v2, " + char(34) + arg + char(34))
	return(st_data((1::nrows), 2))
}

// shortest-path distances (hops) of the symmetrised network; unreachable pairs get the largest finite distance + 1
real matrix lsm_geodesic(real matrix y)
{
	real scalar n, k, mx
	real matrix A, D, T

	n = rows(y)
	A = (editmissing(y, 0) :!= 0)
	A = A :| A'
	D = J(n, n, .)
	D = D :+ 0
	D = (A :* 1) :+ (!A :* 1e9)
	_diag(D, 0)
	for (k = 1; k <= n; k++) {
		T = D[., k] * J(1, n, 1) + J(n, 1, 1) * D[k, .]
		D = D :* (D :<= T) :+ T :* (D :> T)
	}
	mx = max(select(vec(D), vec(D) :< 1e8))
	D = D :* (D :< 1e8) :+ (D :>= 1e8) :* (mx + 1)
	_diag(D, 0)
	return(D)
}

// classical multidimensional scaling: n x d coordinates
real matrix lsm_mds(real matrix D, real scalar d)
{
	real scalar n, k
	real matrix B, J, V
	real rowvector lam

	n = rows(D)
	J = I(n) :- 1 / n
	B = -0.5 * J * (D :* D) * J
	B = (B + B') / 2
	symeigensystem(B, V, lam)
	for (k = 1; k <= d; k++) lam[k] = max((lam[k], 1e-6))
	return(V[., 1..d] :* sqrt(lam[1..d]))
}

// start values: positions (n x d) and coefficients (1 x K) of a logit with the MDS distance as offset
void lsm_start(real matrix y, real scalar dir, real matrix X, real scalar d, real scalar random, real scalar seed,
	real matrix Z0, real rowvector b0, real scalar bil)
{
	real scalar n, nc, k, it, m, i, j
	real matrix Dm, M, XX, S
	real colvector off, yy, bb, p, w, g, keep, r, c, sd

	n = rows(y)
	nc = rows(X) / n
	if (random) {
		S = rseed()
		rseed(seed)
		Z0 = rnormal(n, d, 0, 1)
		rseed(S)
		Z0 = Z0 :- mean(Z0)
	}
	else {
		Dm = lsm_geodesic(y)
		Z0 = lsm_mds(Dm, d)
		if (min(diagonal(variance(Z0))) < 1e-8) {
			S = rseed(); rseed(seed); Z0 = Z0 + rnormal(n, d, 0, 0.1); rseed(S)
		}
	}
	// dyads
	keep = (y :< .)
	_diag(keep, 0)
	if (!dir) keep = uppertriangle(keep, 0)
	m = sum(keep)
	r = selectindex(vec(keep))
	yy = vec(y)[r]
	XX = J(m, 1, 1)
	for (k = 1; k <= nc; k++) XX = XX, vec(X[((k - 1) * n + 1)..(k * n), .])[r]
	Dm = J(n, n, 0)
	for (k = 1; k <= d; k++) Dm = Dm + (Z0[., k] * J(1, n, 1) - J(n, 1, 1) * Z0[., k]') :^ 2
	if (bil) {
		Dm = Z0 * Z0'
		off = vec(Dm)[r]
	}
	else off = -vec(sqrt(Dm))[r]
	bb = J(nc + 1, 1, 0)
	for (it = 1; it <= 30; it++) {
		p = 1 :/ (1 :+ exp(-(XX * bb + off)))
		w = p :* (1 :- p) :+ 1e-9
		g = cross(XX, yy - p) - bb * 1e-3
		M = cross(XX, w, XX) + I(nc + 1) * 1e-3
		bb = bb + cholsolve(M, g)
		bb = rowmax((rowmin((bb, J(nc + 1, 1, 8))), J(nc + 1, 1, -8)))
		if (max(abs(cholsolve(M, g))) < 1e-8) break
	}
	b0 = bb'
}

// start values for a valued network (family 1 Poisson, 2 binomial, 3 normal): positions by classical scaling of a dissimilarity made from the
// symmetrised values (high values = close), coefficients by a few iterations of the GLM of the family with the distance as offset; the start
// of the dispersion (normal family) is the residual variance.
void lsm_vstart(real matrix y, real scalar dir, real matrix X, real scalar d, real scalar random, real scalar seed,
	real matrix Z0, real rowvector b0, real scalar bil, real scalar fam, real matrix tr, real scalar disp0)
{
	real scalar n, nc, k, it, m, mx, mn, p0
	real matrix Dm, M, XX, S, ym
	real colvector off, yy, bb, p, w, g, keep, r, c, tt, eta, mu, res

	n = rows(y)
	nc = rows(X) / n
	if (random) {
		S = rseed()
		rseed(seed)
		Z0 = rnormal(n, d, 0, 1)
		rseed(S)
		Z0 = Z0 :- mean(Z0)
	}
	else {
		ym = editmissing(y, 0)
		if (fam == 2) ym = ym :/ (tr :+ (tr :== 0))
		ym = (ym + ym') / 2
		mx = max(ym)
		mn = min(ym)
		Dm = 1 :+ 2 * (mx :- ym) / (mx - mn + 1e-9)
		_diag(Dm, 0)
		Z0 = lsm_mds(Dm, d)
		if (min(diagonal(variance(Z0))) < 1e-8) {
			S = rseed(); rseed(seed); Z0 = Z0 + rnormal(n, d, 0, 0.1); rseed(S)
		}
	}
	keep = (y :< .)
	_diag(keep, 0)
	if (!dir) keep = uppertriangle(keep, 0)
	m = sum(keep)
	r = selectindex(vec(keep))
	yy = vec(y)[r]
	tt = fam == 2 ? vec(tr)[r] : J(m, 1, 1)
	XX = J(m, 1, 1)
	for (k = 1; k <= nc; k++) XX = XX, vec(X[((k - 1) * n + 1)..(k * n), .])[r]
	Dm = J(n, n, 0)
	for (k = 1; k <= d; k++) Dm = Dm + (Z0[., k] * J(1, n, 1) - J(n, 1, 1) * Z0[., k]') :^ 2
	if (bil) {
		Dm = Z0 * Z0'
		off = vec(Dm)[r]
	}
	else off = -vec(sqrt(Dm))[r]
	bb = J(nc + 1, 1, 0)
	if (fam == 1) bb[1] = ln(max((mean(yy), 0.05)))
	if (fam == 2) {
		p0 = min((max((mean(yy :/ tt), 0.02)), 0.98))
		bb[1] = ln(p0 / (1 - p0))
	}
	for (it = 1; it <= 40; it++) {
		eta = XX * bb + off
		if (fam == 1) {
			eta = rowmax((rowmin((eta, J(m, 1, 12))), J(m, 1, -12)))
			mu = exp(eta)
			w = mu :+ 1e-9
			g = cross(XX, yy - mu) - bb * 1e-3
		}
		else if (fam == 2) {
			p = 1 :/ (1 :+ exp(-eta))
			w = tt :* p :* (1 :- p) :+ 1e-9
			g = cross(XX, yy - tt :* p) - bb * 1e-3
		}
		else {
			w = J(m, 1, 1)
			g = cross(XX, yy - eta) - bb * 1e-3
		}
		M = cross(XX, w, XX) + I(nc + 1) * 1e-3
		bb = bb + cholsolve(M, g)
		bb = rowmax((rowmin((bb, J(nc + 1, 1, 8))), J(nc + 1, 1, -8)))
		if (max(abs(cholsolve(M, g))) < 1e-8) break
	}
	b0 = bb'
	disp0 = 1
	if (fam == 3) {
		res = yy - XX * bb - off
		disp0 = max((variance(res), 1e-3))
	}
}

// checks the values (and trials) of a valued family on the observed dyads; returns 0 if fine, else a code: 1 negative or non-integer value, 2 value above the trials, 3 non-finite value
real scalar lsm_checkvalues(real matrix y, real scalar dir, real scalar fam, real matrix tr)
{
	real scalar n, i, j

	n = rows(y)
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j | (!dir & j < i) | y[i, j] >= .) continue
		if (fam == 3) continue
		if (y[i, j] < 0 | y[i, j] != floor(y[i, j])) return(1)
		if (fam == 2) {
			if (tr[i, j] >= . | tr[i, j] < 1 | tr[i, j] != floor(tr[i, j])) return(4)
			if (y[i, j] > tr[i, j]) return(2)
		}
	}
	return(0)
}

// prior variances: 9 for the intercept; 9 / (mean square of the covariate over the dyads) for a covariate (latentnet's adjust.beta.var)
real rowvector lsm_priorvar(real matrix X, real scalar n, real scalar dir, real scalar basevar)
{
	real scalar nc, k
	real matrix Xk, keep
	real rowvector pv

	nc = rows(X) / n
	pv = J(1, nc + 1, basevar)
	keep = J(n, n, 1)
	_diag(keep, 0)
	if (!dir) keep = uppertriangle(keep, 0)
	for (k = 1; k <= nc; k++) {
		Xk = X[((k - 1) * n + 1)..(k * n), .]
		pv[k + 1] = basevar / (sum(keep :* Xk :^ 2) / sum(keep))
	}
	return(pv)
}

// Monte Carlo standard error of the mean of each column (batch means, batches of about sqrt(S) draws). With C > 1 chains (stored one
// after the other) the standard error of the pooled mean combines the standard errors of the chains.
real rowvector lsm_mcse(real matrix x, | real scalar C)
{
	real scalar S, b, nb, k, c, Sc
	real matrix bm
	real rowvector v

	if (args() < 2) C = 1
	S = rows(x)
	if (C > 1) {
		Sc = S / C
		v = J(1, cols(x), 0)
		for (c = 1; c <= C; c++) v = v + lsm_mcse(x[((c - 1) * Sc + 1)..(c * Sc), .], 1) :^ 2
		return(sqrt(v) / C)
	}
	b = max((1, floor(sqrt(S))))
	nb = floor(S / b)
	if (nb < 2) return(sqrt(diagonal(variance(x))' / S))
	bm = J(nb, cols(x), .)
	for (k = 1; k <= nb; k++) bm[k, .] = mean(x[((k - 1) * b + 1)..(k * b), .])
	return(sqrt(diagonal(variance(bm))' / nb))
}

// Gelman-Rubin potential scale reduction factor of each column; the C chains are stored one after the other
real rowvector lsm_rhat(real matrix x, real scalar C)
{
	real scalar Sc, c, k
	real matrix cm, cv
	real rowvector W, Bn, V

	Sc = rows(x) / C
	cm = J(C, cols(x), .)
	cv = J(C, cols(x), .)
	for (c = 1; c <= C; c++) {
		cm[c, .] = mean(x[((c - 1) * Sc + 1)..(c * Sc), .])
		cv[c, .] = diagonal(variance(x[((c - 1) * Sc + 1)..(c * Sc), .]))'
	}
	W = mean(cv)
	Bn = diagonal(variance(cm))'
	V = (Sc - 1) / Sc :* W + Bn
	return(sqrt(V :/ (W :+ (W :== 0) * 1e-300)))
}

// Procrustes (rotation and reflection, no scaling) of Z onto Ref, both centred: the orthogonal matrix R that makes Z * R closest to Ref
real matrix lsm_procR(real matrix Z, real matrix Ref)
{
	real matrix M, U, Vt
	real rowvector s

	M = cross(Z, Ref)
	svd(M, U, s, Vt)
	return(U * Vt)
}

real matrix lsm_procrustes(real matrix Z, real matrix Ref)
{
	return(Z * lsm_procR(Z, Ref))
}

// all permutations of 1..G, one per row
real matrix lsm_perms(real scalar G)
{
	real matrix P, Q
	real scalar k

	if (G == 1) return(1)
	Q = lsm_perms(G - 1)
	P = J(0, G, .)
	for (k = 1; k <= G; k++) {
		P = P \ ((k > 1 ? Q[., 1..(k - 1)] : J(rows(Q), 0, .)), J(rows(Q), 1, G), (k <= G - 1 ? Q[., k..(G - 1)] : J(rows(Q), 0, .)))
	}
	return(P)
}

// Label correction of the cluster memberships (rows: draws, columns: actors, labels 1..G). The labels of each draw are permuted to agree as
// well as possible with a reference membership matrix, first the one-hot memberships of draw `start`, then (3 rounds) the posterior
// membership probabilities found so far. PI[it, g] = the corrected label of the label g in draw it; Q = n x G posterior membership probabilities.
void lsm_relabel(real matrix Kd, real scalar G, real scalar start, real matrix PI, real matrix Q)
{
	real scalar n, S, rep, it, a, best
	real matrix perms, Qref, M, Qn
	real colvector vM, sc
	real rowvector kr

	n = cols(Kd)
	S = rows(Kd)
	perms = lsm_perms(G)
	Qref = ((Kd[start, .]' * J(1, G, 1)) :== (J(n, 1, 1) * (1..G))) :* 1
	PI = J(S, G, .)
	Q = Qref
	for (rep = 1; rep <= 3; rep++) {
		Qn = J(n, G, 0)
		for (it = 1; it <= S; it++) {
			M = ((Kd[it, .]' * J(1, G, 1)) :== (J(n, 1, 1) * (1..G)))' * Qref
			vM = vec(M)
			sc = J(rows(perms), 1, 0)
			for (a = 1; a <= G; a++) sc = sc + vM[(perms[., a] :- 1) * G :+ a]
			best = selectindex(sc :== max(sc))[1]
			PI[it, .] = perms[best, .]
			kr = PI[it, Kd[it, .]]
			Qn = Qn + ((kr' * J(1, G, 1)) :== (J(n, 1, 1) * (1..G)))
		}
		Qref = Qn / S
		Q = Qref
	}
}



// n x n tie matrix of a network object's matrix: 0/1, missing stays missing, self-loops are ignored later
real matrix lsm_ymat(real matrix a)
{
	real matrix y

	y = (editmissing(a, 0) :!= 0) :/ (a :< .)
	return(y)
}

// one numeric variable per actor -> |x_i - x_j| matrix
real matrix lsm_node_covs(string scalar var, real scalar n)
{
	real colvector v

	v = st_data(., var)
	return(abs(v * J(1, n, 1) - J(n, 1, 1) * v'))
}

// effective sample size of a series (Geyer's initial monotone positive sequence)
real scalar lsm_ess(real colvector x)
{
	real scalar S, m, k, g0, gam, prev, sumg, G
	real colvector xc

	S = rows(x)
	xc = x :- mean(x)
	g0 = cross(xc, xc) / S
	if (g0 <= 0) return(S)
	sumg = 0
	prev = .
	for (k = 0; k < S - 2; k = k + 2) {
		G = (cross(xc[1..(S - k)], xc[(k + 1)..S]) + cross(xc[1..(S - k - 1)], xc[(k + 2)..S])) / S
		if (G <= 0) break
		if (prev < . & G > prev) G = prev
		prev = G
		sumg = sumg + G
	}
	// sum over pairs from lag 0: tau = -g0 + 2 sum_m (g_2m + g_2m+1)
	gam = (-g0 + 2 * sumg) / g0
	if (gam < 1 / S) gam = 1 / S
	return(min((S, S / gam)))
}

// Geweke z-score: first 10% of the series against the last 50%; the variance of each mean uses the effective sample size of its segment
real scalar lsm_geweke(real colvector x)
{
	real scalar S, na, nb, va, vb
	real colvector a, b

	S = rows(x)
	na = max((2, floor(0.1 * S)))
	nb = max((2, floor(0.5 * S)))
	a = x[1..na]
	b = x[(S - nb + 1)..S]
	va = variance(a) / lsm_ess(a)
	vb = variance(b) / lsm_ess(b)
	if (va + vb <= 0) return(0)
	return((mean(a) - mean(b)) / sqrt(va + vb))
}

// diagnostics of the columns of a draws matrix (C chains stored one after the other): (ESS, Geweke z, mean, sd, Gelman-Rubin R) per column.
// ESS: sum over the chains; Geweke z: the value of the largest absolute size over the chains; R: missing for one chain.
real matrix LsmDiag(real matrix x, | real scalar C)
{
	real scalar k, c, Sc, z
	real matrix R
	real colvector xc

	if (args() < 2) C = 1
	Sc = rows(x) / C
	R = J(cols(x), 5, .)
	for (k = 1; k <= cols(x); k++) {
		R[k, 1] = 0
		R[k, 2] = 0
		for (c = 1; c <= C; c++) {
			xc = x[((c - 1) * Sc + 1)..(c * Sc), k]
			R[k, 1] = R[k, 1] + lsm_ess(xc)
			z = lsm_geweke(xc)
			if (abs(z) > abs(R[k, 2])) R[k, 2] = z
		}
		R[k, 3] = mean(x[., k])
		R[k, 4] = sqrt(variance(x[., k]))
		if (C > 1) R[k, 5] = lsm_rhat(x[., k], C)
	}
	return(R)
}

// summary of a variance parameter: mean, sd, mcse, 2.5%, 97.5%
real rowvector lsm_vartable(real colvector v, | real scalar C)
{
	real colvector srt
	real scalar S

	if (args() < 2) C = 1
	S = rows(v)
	srt = sort(v, 1)
	return((mean(v), sqrt(variance(v)), lsm_mcse(v, C), srt[max((1, ceil(0.025 * S)))], srt[min((S, floor(0.975 * S) + 1))]))
}

// centre the positions and rotate them to their principal axes (a rigid motion; positions are identified up to such motions only)
real matrix lsm_ztrans(real matrix Z, real scalar bil)
{
	real matrix Zc, V, U, Vt
	real rowvector sv

	// the bilinear model is not invariant to translations: positions are only rotated or reflected there
	Zc = bil ? Z : Z :- mean(Z)
	if (cols(Zc) > 1) {
		svd(Zc, U, sv, Vt)
		Zc = Zc * Vt'
	}
	return(Zc)
}

/*
	LsmFit: runs the chain and stores the summary. The posterior sample stays in the Mata externals __nwlsm_beta (S x K), __nwlsm_llk,
	__nwlsm_s2 (S x 1), __nwlsm_ev (S x 3: variances of the sender, receiver and sociality effects), __nwlsm_Z (S x n d, by actor then
	dimension), __nwlsm_E (S x 3 n: the effects, by type then actor), __nwlsm_D (mean latent term) and __nwlsm_P (mean tie probabilities).
	opt = (sender, receiver, sociality, bilinear, var_s, df_s, var_r, df_r, var_o, df_o, mklstarts); nm = the names of the Stata matrices to fill:
	B V T D MKLZ A P MB ZM ME VT ES ER EO (means of the draws, covariance, coefficient table, mean distances, MKL positions, tuning,
	variance table, effects (n x 2: mean, sd), mean tie probabilities, MKL coefficients, mean of aligned draws, MKL effects).
	The order of nm: B V T D MKLZ ACC P MB ZM ME VT ES ER EO.
*/
void LsmFit(real matrix y, real scalar dir, real matrix X, real scalar d, real scalar nsample, real scalar burnin,
	real scalar interval, real scalar seed, real scalar pilot, real scalar plen, real scalar random,
	real rowvector opt, string rowvector nm, | real matrix tr, real matrix io)
{
	class LsmNat scalar N
	real scalar fam, disp0, n, nc, K, S, Sc, C, G, ng, k, it, rows_, off, nd, matok, na, t, i, g2, bi, hasio, q, s20
	real rowvector prior, pvar, b0, flags, mb, act, cl, isv
	real matrix Ob, Wm, PR, Bj
	real colvector Os2, Ol, bols, ores
	real matrix Z0, DispD, beta, Zs, Dbar, Zm, Za, Zr, T, EP, Pm, VT, Es, Ev, Eall, EPrior, E0, Em, S2, RH, Mu, Lam, Kd, PI, Q, S2r, Lr, Kr, MUa, R, Zi
	real colvector res, llk, srt, shift, cls
	string scalar Dn
	string rowvector vn

	n = rows(y)
	nc = rows(X) / n
	K = nc + 1
	C = opt[12]
	G = opt[14]
	ng = max((G, 1))
	Sc = ceil(nsample / C)
	S = Sc * C
	nd = n * d
	flags = opt[1..4]
	act = opt[1..3]
	na = sum(act)
	pvar = lsm_priorvar(X, n, dir, 9)
	prior = (0.125 * n^(2 / d), sqrt(n), 0.125 * n^(2 / d))
	cl = J(1, 5, 0)
	if (G > 0) cl = (G, 0.125 * (n / G)^(2 / d), sqrt(n / G), 2 * 0.125 * n^(2 / d), sqrt(n / G))
	EPrior = (opt[5], opt[6], opt[5] \ opt[7], opt[8], opt[7] \ opt[9], opt[10], opt[9])
	E0 = J(n, 3, 0)
	if (args() < 14) tr = J(0, 0, .)
	fam = (args() >= 14 ? opt[15] : 0)
	disp0 = 1
	if (fam == 0) lsm_start(y, dir, X, d, random, seed, Z0, b0, flags[4])
	else lsm_vstart(y, dir, X, d, random, seed, Z0, b0, flags[4], fam, tr, disp0)
	if (args() < 15) io = J(0, 0, .)
	hasio = (rows(io) > 0)
	q = 0
	isv = J(1, 6, 0)
	if (hasio) {
		// start: least squares of the outcome on the regressors without influence; influence effect 0
		q = cols(io) - 1
		bols = invsym(cross(io[., 3..cols(io)], io[., 3..cols(io)])) * cross(io[., 3..cols(io)], io[., 1])
		ores = io[., 1] - io[., 3..cols(io)] * bols
		s20 = max((cross(ores, ores) / max((n - q + 1, 1)), 1e-4))
		isv = (opt[18], opt[19], opt[20], opt[21], opt[22], s20, bols', 0)
	}
	N.open()
	N.setup(y, dir, X, d, prior, pvar, b0, Z0, flags, EPrior, E0, cl, (fam, (fam == 3 ? opt[16] : 0), (fam == 3 ? opt[17] : 0), (fam == 3 ? disp0 : 0)), tr, io, isv)
	rows_ = 2 * K + 8 + S * K + S + S + S * ng + 3 * S + S * nd + na * S * n + (G > 0 ? S * G * d + S * G + S * n : 0) + 2 * n * n + nd + K + na * n + 3 + (hasio ? S * q + 2 * S + n * n + 3 * n : 0)
	res = N.call(sprintf("10 %f %f %f %f %f %f %f %f %f %f %f %f %f %f", N.L, seed, burnin, interval, Sc, pilot, plen, 0.6, 0.6, 0.234, 0.6, opt[11], C, opt[13]), rows_)
	N.close()
	if (hasmissing(res)) {
		errprintf("nwlsm: the sampler returned missing values\n")
		exit(498)
	}
	off = 2 * K + 8
	beta = colshape(res[(off + 1)..(off + S * K)], K)
	off = off + S * K
	llk = res[(off + 1)..(off + S)]
	off = off + S
	DispD = res[(off + 1)..(off + S)]
	off = off + S
	S2 = colshape(res[(off + 1)..(off + S * ng)], ng)
	off = off + S * ng
	Ev = colshape(res[(off + 1)..(off + 3 * S)]', S)'
	off = off + 3 * S
	Zs = colshape(res[(off + 1)..(off + S * nd)], nd)
	off = off + S * nd
	Eall = J(S, 3 * n, 0)
	for (t = 1; t <= 3; t++) {
		if (act[t]) {
			Eall[., ((t - 1) * n + 1)..(t * n)] = colshape(res[(off + 1)..(off + S * n)], n)
			off = off + S * n
		}
	}
	if (G > 0) {
		Mu = colshape(res[(off + 1)..(off + S * G * d)], G * d)
		off = off + S * G * d
		Lam = colshape(res[(off + 1)..(off + S * G)], G)
		off = off + S * G
		Kd = colshape(res[(off + 1)..(off + S * n)], n)
		off = off + S * n
	}
	Dbar = colshape(res[(off + 1)..(off + n * n)], n)
	off = off + n * n
	Pm = colshape(res[(off + 1)..(off + n * n)], n)
	off = off + n * n
	Zm = colshape(res[(off + 1)..(off + nd)]', d)
	off = off + nd
	mb = res[(off + 1)..(off + K)]'
	off = off + K
	Em = J(n, 3, 0)
	for (t = 1; t <= 3; t++) {
		if (act[t]) {
			Em[., t] = res[(off + 1)..(off + n)]
			off = off + n
		}
	}
	st_local("__mklstart", strofreal(res[off + 1]))
	st_local("__mkliter", strofreal(res[off + 2]))
	st_local("__mklobj", strofreal(res[off + 3], "%18.10g"))
	bi = res[off + 1]
	off = off + 3
	if (hasio) {
		Ob = colshape(res[(off + 1)..(off + S * q)], q)
		off = off + S * q
		Os2 = res[(off + 1)..(off + S)]
		off = off + S
		Ol = res[(off + 1)..(off + S)]
		off = off + S
		Wm = colshape(res[(off + 1)..(off + n * n)], n)
		off = off + n * n
		PR = io[., 1], res[(off + 1)..(off + n)], res[(off + n + 1)..(off + 2 * n)], res[(off + 2 * n + 1)..(off + 3 * n)]
		off = off + 3 * n
	}
	// MKL: intercept absorbs the means of the effects (sociality enters twice), the positions are centred
	shift = 0
	if (act[1]) {
		shift = shift + mean(Em[., 1])
		Em[., 1] = Em[., 1] :- mean(Em[., 1])
	}
	if (act[2]) {
		shift = shift + mean(Em[., 2])
		Em[., 2] = Em[., 2] :- mean(Em[., 2])
	}
	if (act[3]) {
		shift = shift + 2 * mean(Em[., 3])
		Em[., 3] = Em[., 3] :- mean(Em[., 3])
	}
	mb[1] = mb[1] + shift
	Zm = lsm_ztrans(Zm, flags[4])
	// label correction of the clusters
	if (G > 0) {
		lsm_relabel(Kd, G, bi, PI, Q)
		S2r = J(S, G, .)
		Lr = J(S, G, .)
		for (it = 1; it <= S; it++) {
			for (g2 = 1; g2 <= G; g2++) {
				S2r[it, PI[it, g2]] = S2[it, g2]
				Lr[it, PI[it, g2]] = Lam[it, g2]
			}
			Kd[it, .] = PI[it, Kd[it, .]]
		}
		S2 = S2r
		Lam = Lr
	}
	// summary of the coefficients
	Bj = hasio ? (beta, Ob) : beta
	T = J(K + q, 5, .)
	T[., 1] = mean(Bj)'
	T[., 2] = sqrt(diagonal(variance(Bj)))
	T[., 3] = lsm_mcse(Bj, C)'
	for (k = 1; k <= K + q; k++) {
		srt = sort(Bj[., k], 1)
		T[k, 4] = srt[max((1, ceil(0.025 * S)))]
		T[k, 5] = srt[min((S, floor(0.975 * S) + 1))]
	}
	st_matrix(nm[1], mean(Bj))
	st_matrix(nm[2], variance(Bj))
	st_matrix(nm[3], T)
	// variance parameters: sigma2 of the positions (one per cluster with clusters), then the active effects
	VT = J(0, 5, .)
	for (k = 1; k <= ng; k++) VT = VT \ lsm_vartable(S2[., k], C)
	for (t = 1; t <= 3; t++) if (act[t]) VT = VT \ lsm_vartable(Ev[., t], C)
	if (fam == 3) VT = VT \ lsm_vartable(DispD, C)
	if (hasio) VT = VT \ lsm_vartable(Os2, C)
	st_matrix(nm[11], VT)
	// Gelman-Rubin statistics (several chains): coefficients, log likelihood, variance parameters
	if (C > 1) {
		RH = lsm_rhat(Bj, C)' \ lsm_rhat(llk, C)
		for (k = 1; k <= ng; k++) RH = RH \ lsm_rhat(S2[., k], C)
		for (t = 1; t <= 3; t++) if (act[t]) RH = RH \ lsm_rhat(Ev[., t], C)
		if (fam == 3) RH = RH \ lsm_rhat(DispD, C)
		if (hasio) RH = RH \ lsm_rhat(Os2, C)
		st_matrix(nm[15], RH)
	}
	// positions: the mean of the draws, each aligned (rotation, reflection) to the MKL positions; the cluster means follow the same rotations
	Zr = Zm
	Za = J(n, d, 0)
	MUa = J(G, d, 0)
	for (it = 1; it <= S; it++) {
		Zi = colshape(Zs[it, .], d)
		if (!flags[4]) {
			if (G > 0) {
				Mu[it, .] = vec((colshape(Mu[it, .], d) :- mean(Zi))')'
			}
			Zi = Zi :- mean(Zi)
		}
		R = lsm_procR(Zi, Zr)
		Za = Za + Zi * R
		if (G > 0) {
			Zi = colshape(Mu[it, .], d) * R
			for (g2 = 1; g2 <= G; g2++) MUa[PI[it, g2], .] = MUa[PI[it, g2], .] + Zi[g2, .] / S
		}
	}
	Za = Za / S
	// effects: posterior mean and sd
	matok = (n <= c("max_matdim"))
	for (t = 1; t <= 3; t++) {
		if (act[t]) {
			Es = (mean(Eall[., ((t - 1) * n + 1)..(t * n)])', sqrt(diagonal(variance(Eall[., ((t - 1) * n + 1)..(t * n)]))))
			if (matok) st_matrix(nm[11 + t], Es)
		}
	}
	if (matok) {
		st_matrix(nm[4], Dbar)
		st_matrix(nm[5], Zm)
		st_matrix(nm[7], Pm)
		st_matrix(nm[9], Za)
		st_matrix(nm[10], Em)
		if (G > 0) st_matrix(nm[16], Q)
		if (hasio) {
			st_matrix(nm[20], Wm)
			st_matrix(nm[21], PR)
		}
	}
	if (G > 0) {
		// hard classification: the most probable cluster (ties: the lower number)
		st_matrix(nm[17], MUa)
		st_matrix(nm[18], (mean(Lam)', sqrt(diagonal(variance(Lam)))))
		if (matok) {
			cls = J(n, 1, .)
			for (i = 1; i <= n; i++) cls[i] = selectindex(Q[i, .] :== max(Q[i, .]))[1]
			st_matrix(nm[19], cls)
		}
	}
	st_matrix(nm[8], mb)
	st_matrix(nm[6], (res[1..(K + 1)]', res[(K + 2)..(K + 4)]' \ res[(K + 5)..(2 * K + 5)]', res[(2 * K + 6)..(2 * K + 8)]'))
	st_local("__matok", strofreal(matok))
	st_local("__llk", strofreal(mean(llk), "%18.10g"))
	st_local("__llksd", strofreal(sqrt(variance(llk)), "%18.10g"))
	st_local("__s2", strofreal(mean(S2[., 1]), "%18.10g"))
	st_local("__s2sd", strofreal(sqrt(variance(S2[., 1])), "%18.10g"))
	st_local("__nsample_total", strofreal(S))
	if (hasio) {
		st_local("__pb2", strofreal(mean(Ob[., q] :> 0), "%18.10g"))
		st_local("__ollk", strofreal(mean(Ol), "%18.10g"))
		st_local("__ollksd", strofreal(sqrt(variance(Ol)), "%18.10g"))
	}
	st_local("__priorscale", strofreal(G > 0 ? cl[2] : prior[1], "%18.10g"))
	st_local("__priordf", strofreal(G > 0 ? cl[3] : prior[2], "%18.10g"))
	st_matrix(nm[4] + "_pv", pvar)
	st_matrix(nm[4] + "_b0", b0)
	external __nwlsm_beta, __nwlsm_llk, __nwlsm_s2, __nwlsm_Z, __nwlsm_D, __nwlsm_ZA, __nwlsm_ev, __nwlsm_E, __nwlsm_P, __nwlsm_MZ, __nwlsm_C, __nwlsm_K, __nwlsm_lam, __nwlsm_disp, __nwlsm_fam
	__nwlsm_disp = DispD
	__nwlsm_fam = fam
	external __nwlsm_ob, __nwlsm_os2, __nwlsm_ol, __nwlsm_io
	__nwlsm_ob = (hasio ? Ob : J(S, 0, .))
	__nwlsm_os2 = (hasio ? Os2 : J(S, 0, .))
	__nwlsm_ol = (hasio ? Ol : J(0, 0, .))
	__nwlsm_io = io
	external __nwlsm_gY, __nwlsm_gX, __nwlsm_gT, __nwlsm_gm
	__nwlsm_gY = y
	__nwlsm_gX = X
	__nwlsm_gT = tr
	__nwlsm_gm = (dir, d, flags[1], flags[2], flags[3], flags[4], fam)
	__nwlsm_beta = beta
	__nwlsm_llk = llk
	__nwlsm_s2 = S2
	__nwlsm_Z = Zs
	__nwlsm_D = Dbar
	__nwlsm_ZA = Za
	__nwlsm_ev = Ev
	__nwlsm_E = Eall
	__nwlsm_P = Pm
	__nwlsm_MZ = Zm
	__nwlsm_C = C
	__nwlsm_K = (G > 0 ? Kd : J(0, 0, .))
	__nwlsm_lam = (G > 0 ? Lam : J(0, 0, .))
}

// the posterior sample kept in Mata by the last fit: LsmSample("beta"), "llk", "s2", "ev", "E", "Z" (one row per draw), "D", "P", "Zalign", "mkl"
transmorphic matrix LsmSample(string scalar what)
{
	external __nwlsm_beta, __nwlsm_llk, __nwlsm_s2, __nwlsm_Z, __nwlsm_D, __nwlsm_ZA, __nwlsm_ev, __nwlsm_E, __nwlsm_P, __nwlsm_MZ, __nwlsm_C, __nwlsm_K, __nwlsm_lam, __nwlsm_disp, __nwlsm_fam

	external __nwlsm_ob, __nwlsm_os2, __nwlsm_ol, __nwlsm_io
	if (what == "ob") return(__nwlsm_ob)
	if (what == "os2") return(__nwlsm_os2)
	if (what == "ollk") return(__nwlsm_ol)
	if (what == "io") return(__nwlsm_io)
	if (what == "disp") return(__nwlsm_disp)
	if (what == "family") return(__nwlsm_fam)
	if (what == "chains") return(__nwlsm_C)
	if (what == "class") return(__nwlsm_K)
	if (what == "lambda") return(__nwlsm_lam)
	if (what == "beta") return(__nwlsm_beta)
	if (what == "llk") return(__nwlsm_llk)
	if (what == "s2") return(__nwlsm_s2)
	if (what == "ev") return(__nwlsm_ev)
	if (what == "E") return(__nwlsm_E)
	if (what == "Z") return(__nwlsm_Z)
	if (what == "D") return(__nwlsm_D)
	if (what == "P") return(__nwlsm_P)
	if (what == "Zalign") return(__nwlsm_ZA)
	if (what == "mkl") return(__nwlsm_MZ)
	errprintf("LsmSample(): unknown item\n")
	exit(198)
}

// writes the draws of the coefficients, the log likelihood and the variances into the (empty) current dataset
void LsmSaveDraws(string rowvector names, real rowvector act, | string rowvector onames)
{
	real scalar k, S
	string scalar vn
	string rowvector vnn
	external __nwlsm_beta, __nwlsm_llk, __nwlsm_s2, __nwlsm_ev

	S = rows(__nwlsm_beta)
	st_addobs(S)
	for (k = 1; k <= cols(__nwlsm_beta); k++) {
		vn = (names[k] == "_cons" ? "beta_cons" : names[k])
		(void) st_addvar("double", vn)
		st_store(., vn, __nwlsm_beta[., k])
	}
	(void) st_addvar("double", "loglik")
	st_store(., "loglik", __nwlsm_llk)
	if (cols(__nwlsm_s2) == 1) {
		(void) st_addvar("double", "sigma2")
		st_store(., "sigma2", __nwlsm_s2)
	}
	else {
		for (k = 1; k <= cols(__nwlsm_s2); k++) {
			(void) st_addvar("double", "sigma2_" + strofreal(k))
			st_store(., "sigma2_" + strofreal(k), __nwlsm_s2[., k])
		}
		external __nwlsm_lam
		for (k = 1; k <= cols(__nwlsm_lam); k++) {
			(void) st_addvar("double", "lambda_" + strofreal(k))
			st_store(., "lambda_" + strofreal(k), __nwlsm_lam[., k])
		}
	}
	vnn = ("sender_var", "receiver_var", "sociality_var")
	for (k = 1; k <= 3; k++) {
		if (act[k]) {
			(void) st_addvar("double", vnn[k])
			st_store(., vnn[k], __nwlsm_ev[., k])
		}
	}
	external __nwlsm_disp, __nwlsm_fam
	if (__nwlsm_fam == 3) {
		(void) st_addvar("double", "dispersion")
		st_store(., "dispersion", __nwlsm_disp)
	}
	if (args() >= 3 & cols(LsmSample("ob")) > 0) {
		for (k = 1; k <= cols(onames); k++) {
			(void) st_addvar("double", "out_" + onames[k])
			st_store(., "out_" + onames[k], LsmSample("ob")[., k])
		}
		(void) st_addvar("double", "sigma2_y")
		st_store(., "sigma2_y", LsmSample("os2"))
		(void) st_addvar("double", "loglik_outcome")
		st_store(., "loglik_outcome", LsmSample("ollk"))
	}
}

// summary matrix of one statistic: obs, min, mean, max, MC p-value (as latentnet: 2 min(P(sim <= obs), P(sim >= obs)), at most 1); rows: cells
real matrix lsm_gofsum(real rowvector obs, real matrix sim)
{
	real matrix R
	real scalar k, ns

	ns = rows(sim)
	R = J(cols(obs), 5, .)
	for (k = 1; k <= cols(obs); k++) {
		R[k, 1] = obs[k]
		R[k, 2] = min(sim[., k])
		R[k, 3] = mean(sim[., k])
		R[k, 4] = max(sim[., k])
		R[k, 5] = min((1, 2 * min((mean(sim[., k] :<= obs[k]), mean(sim[., k] :>= obs[k])))))
	}
	return(R)
}

// Goodness of fit: nsim networks simulated from posterior draws (evenly spaced over the stored draws), the statistics of the groups in mask
// (bits: 1 degree or in-degree, 2 out-degree, 4 edgewise shared partners, 8 dyadwise shared partners, 16 geodesic distances, 32 model statistics,
// 64 tie values) for the observed and the simulated networks. Fills the Stata matrices <nm>_<group> (cells x 5: observed, min, mean, max, MC p-value)
// and <nm>_q_<group> (cells x 6: observed share, min, 25%, median, 75%, max of the share over the simulations).
void LsmGof(real scalar nsim, real scalar seed, real scalar mask, real scalar threads, real scalar M, string scalar nm, | real scalar usenet, real scalar thr)
{
	class LsmNat scalar N
	real scalar n, nc, K, d, dir, fam, S, na, k, t, RW, m, ncut, j, off, w, g
	real matrix Dr, E, out, Sm, Q, X, y, tr
	real rowvector flags, act, w7, idx, cuts, ob
	real colvector yv, res, srt
	external __nwlsm_beta, __nwlsm_Z, __nwlsm_E, __nwlsm_disp, __nwlsm_gY, __nwlsm_gX, __nwlsm_gT, __nwlsm_gm
	string rowvector grp

	y = __nwlsm_gY
	X = __nwlsm_gX
	tr = __nwlsm_gT
	n = rows(y)
	nc = rows(X) / n
	K = nc + 1
	dir = __nwlsm_gm[1]
	d = __nwlsm_gm[2]
	flags = __nwlsm_gm[3..6]
	fam = __nwlsm_gm[7]
	act = flags[1..3]
	na = sum(act)
	S = rows(__nwlsm_beta)
	if (nsim > S) nsim = S
	idx = (nsim == 1 ? 1 : round(1 :+ (0..(nsim - 1)) * (S - 1) / (nsim - 1)))
	Dr = __nwlsm_beta[idx', .], __nwlsm_Z[idx', .]
	for (t = 1; t <= 3; t++) if (act[t]) Dr = Dr, __nwlsm_E[idx', ((t - 1) * n + 1)..(t * n)]
	Dr = Dr, (fam == 3 ? __nwlsm_disp[idx'] : J(nsim, 1, 0))
	cuts = J(1, 0, .)
	ncut = 0
	if (fam == 3 & mod(floor(mask / 64), 2) == 1) {
		// off-diagonal values of the observed network (undirected: each tie once)
		yv = J(0, 1, .)
		for (j = 1; j <= n; j++) {
			if (dir) yv = yv \ select(y[j, .]', (y[j, .]' :< .) :& ((1::n) :!= j))
			else if (j < n) yv = yv \ select(y[j, (j + 1)..n]', y[j, (j + 1)..n]' :< .)
		}
		srt = sort(yv, 1)
		cuts = J(1, 9, .)
		for (k = 1; k <= 9; k++) cuts[k] = srt[max((1, ceil(k * rows(srt) / 10)))]
		ncut = 9
	}
	N.open()
	N.setup(y, dir, X, d, (1, 1, 1), J(1, K, 9), J(1, K, 0), J(n, d, 0), flags, (1, 1, 1 \ 1, 1, 1 \ 1, 1, 1), J(n, 3, 0), J(1, 5, 0), (fam, 1, 1, (fam == 3 ? 1 : 0)), tr)
	N.append(cuts' \ vec(Dr'))
	w7 = J(1, 7, 0)
	if (mod(floor(mask / 1), 2) == 1) w7[1] = n
	if (mod(floor(mask / 2), 2) == 1) w7[2] = n
	if (mod(floor(mask / 4), 2) == 1) w7[3] = n - 1
	if (mod(floor(mask / 8), 2) == 1) w7[4] = n - 1
	if (mod(floor(mask / 16), 2) == 1) w7[5] = n
	if (mod(floor(mask / 32), 2) == 1) w7[6] = 1 + nc
	if (mod(floor(mask / 64), 2) == 1) w7[7] = (fam == 3 ? ncut + 1 : M + 1)
	RW = sum(w7)
	if (args() < 8) {
		usenet = (fam != 3)
		thr = 0
	}
	res = N.call(sprintf("20 %f %f %f %f %f %f %f %f %s", N.L, nsim, seed, threads, mask, M, ncut, usenet, strofreal(thr, "%20.15g")), (nsim + 1) * RW)
	N.close()
	if (hasmissing(res)) {
		errprintf("nwlsm: the simulation returned missing values\n")
		exit(498)
	}
	out = colshape(res, RW)
	grp = ("deg", "odeg", "esp", "dsp", "dist", "model", "val")
	if (dir) grp[1] = "ideg"
	off = 0
	for (g = 1; g <= 7; g++) {
		w = w7[g]
		if (w == 0) continue
		ob = out[1, (off + 1)..(off + w)]
		Sm = out[2..(nsim + 1), (off + 1)..(off + w)]
		st_matrix(nm + "_" + grp[g], lsm_gofsum(ob, Sm))
		if (g != 6) {
			Q = Sm :/ (rowsum(Sm) :+ (rowsum(Sm) :== 0))
			ob = ob / max((sum(ob), 1))
			Q = Q'
			Q = ob', rowmin(Q), J(w, 1, .), J(w, 1, .), J(w, 1, .), rowmax(Q)
			for (k = 1; k <= w; k++) {
				srt = sort(Sm[., k] :/ (rowsum(Sm) :+ (rowsum(Sm) :== 0)), 1)
				Q[k, 3] = srt[max((1, ceil(0.25 * nsim)))]
				Q[k, 4] = srt[max((1, ceil(0.5 * nsim)))]
				Q[k, 5] = srt[max((1, ceil(0.75 * nsim)))]
			}
			st_matrix(nm + "_q_" + grp[g], Q)
		}
		off = off + w
	}
	st_numscalar(nm + "_nsim", nsim)
}

// Posterior predictive check of the expected values: observed dyads (undirected: each once) sorted by the posterior mean prediction and cut into
// G groups of about equal size: (dyads, mean predicted, mean observed) per group in tn; summary in sn: (dyads, Brier score or RMSE, AUC or correlation,
// mean log loss or mean absolute error).
void LsmPPCheck(real scalar G, string scalar tn, string scalar sn)
{
	external __nwlsm_P, __nwlsm_gY, __nwlsm_gm
	real matrix Y, P
	real colvector yv, pv, ord, r1, r0
	real scalar n, dir, fam, i, j, g, lo, hi, m, ng, auc, n1, n0, t, br
	real matrix T

	Y = __nwlsm_gY
	P = __nwlsm_P
	n = rows(Y)
	dir = __nwlsm_gm[1]
	fam = __nwlsm_gm[7]
	yv = J(0, 1, .)
	pv = J(0, 1, .)
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j | (!dir & j < i) | Y[i, j] >= .) continue
		yv = yv \ Y[i, j]
		pv = pv \ P[i, j]
	}
	m = rows(yv)
	ord = order(pv, 1)
	yv = yv[ord]
	pv = pv[ord]
	ng = min((G, m))
	T = J(ng, 3, .)
	for (g = 1; g <= ng; g++) {
		lo = floor((g - 1) * m / ng) + 1
		hi = floor(g * m / ng)
		T[g, .] = (hi - lo + 1, mean(pv[lo..hi]), mean(yv[lo..hi]))
	}
	st_matrix(tn, T)
	if (fam == 0) {
		br = mean((yv - pv) :^ 2)
		// AUC by ranks (ties get mid ranks)
		r1 = J(m, 1, .)
		i = 1
		while (i <= m) {
			j = i
			while (j < m) {
				if (pv[j + 1] != pv[i]) break
				j++
			}
			r1[i..j] = J(j - i + 1, 1, (i + j) / 2)
			i = j + 1
		}
		n1 = sum(yv :== 1)
		n0 = m - n1
		auc = .
		if (n1 > 0 & n0 > 0) auc = (sum(select(r1, yv :== 1)) - n1 * (n1 + 1) / 2) / (n1 * n0)
		st_matrix(sn, (m, br, auc, -mean(yv :* ln(pv :+ (pv :== 0) * 1e-300) + (1 :- yv) :* ln(1 :- pv :+ (pv :== 1) * 1e-300))))
	}
	else {
		auc = .
		if (variance(yv)[1, 1] > 0 & variance(pv)[1, 1] > 0) auc = correlation((yv, pv))[1, 2]
		st_matrix(sn, (m, sqrt(mean((yv - pv) :^ 2)), auc, mean(abs(yv - pv))))
	}
}

// the draws of the coefficients, the log likelihood and the variance parameters (columns in the order of the rows of estat mcmcdiag), chains one after the other
real matrix LsmTraceMat(real rowvector act)
{
	external __nwlsm_beta, __nwlsm_llk, __nwlsm_s2, __nwlsm_ev, __nwlsm_disp, __nwlsm_fam, __nwlsm_ob, __nwlsm_os2
	real matrix M
	real scalar t

	M = (__nwlsm_beta, __nwlsm_ob, __nwlsm_llk, __nwlsm_s2)
	for (t = 1; t <= 3; t++) if (act[t]) M = (M, __nwlsm_ev[., t])
	if (__nwlsm_fam == 3) M = (M, __nwlsm_disp)
	if (cols(__nwlsm_os2) > 0) M = (M, __nwlsm_os2)
	return(M)
}

// test helper: the plugin's log likelihood and log prior of Z at (beta, Z, s2), and its incremental log-likelihood change
// if actor i moves to znew. Returns (loglik, logprior_Z, dll, dll_effect), the last for the effect of type t (0 = none) of actor i
// changing by delta. E: n x 3 start values of the effects; opt as for LsmFit.
real rowvector LsmCheck(real matrix y, real scalar dir, real matrix X, real scalar d, real rowvector beta, real matrix Z,
	real scalar s2, real scalar i, real rowvector znew, real rowvector opt, real matrix E, real scalar t, real scalar delta, | real rowvector fm, real matrix tr)
{
	class LsmNat scalar N
	real colvector r1, r2, r3
	real rowvector prior, pvar
	real matrix EPrior
	real scalar n

	n = rows(y)
	pvar = lsm_priorvar(X, n, dir, 9)
	prior = (1, 1, s2)
	EPrior = (opt[5], opt[6], opt[5] \ opt[7], opt[8], opt[7] \ opt[9], opt[10], opt[9])
	N.open()
	if (args() < 14) fm = (0, 0, 0, 0)
	if (args() < 15) tr = J(0, 0, .)
	N.setup(y, dir, X, d, prior, pvar, beta, Z, opt[1..4], EPrior, E, J(1, 5, 0), fm, tr)
	r1 = N.call(sprintf("1 %f", N.L), 2)
	r2 = N.call(sprintf("2 %f %f ", N.L, i) + invtokens(strofreal(znew, "%20.15g")), 1)
	r3 = 0
	if (t > 0) r3 = N.call(sprintf("3 %f %f %f %s", N.L, t - 1, i, strofreal(delta, "%20.15g")), 1)
	N.close()
	return((r1', r2, r3))
}

// test helper (influence model): outcome log likelihood at the state, the weighted means m_k, then the change of the outcome log likelihood
// when actor i moves to znew by the incremental update, by full recomputation of the moved state, and by the plugin's state after the move
real rowvector LsmCheckInfl(real matrix y, real scalar dir, real matrix X, real scalar d, real rowvector beta, real matrix Z, real scalar s2,
	real matrix io, real rowvector isv, real scalar i, real rowvector znew)
{
	class LsmNat scalar N
	real colvector r1, r2
	real rowvector pvar
	real matrix EPrior
	real scalar n

	n = rows(y)
	pvar = lsm_priorvar(X, n, dir, 9)
	EPrior = (1, 3, 1 \ 1, 3, 1 \ 1, 3, 1)
	N.open()
	N.setup(y, dir, X, d, (1, 1, s2), pvar, beta, Z, (0, 0, 0, 0), EPrior, J(n, 3, 0), J(1, 5, 0), (0, 0, 0, 0), J(0, 0, .), io, isv)
	r1 = N.call(sprintf("5 %f", N.L), n + 1)
	r2 = N.call(sprintf("6 %f %f ", N.L, i) + invtokens(strofreal(znew, "%20.15g")), 3)
	N.close()
	return((r1', r2'))
}

// test helper: objective of the MKL optimisation at the state (target = 0.9 for a tie, 0.15 for none), the analytic gradient and a
// central-difference gradient. Returns (f, analytic' , numeric') as a row vector.
real rowvector LsmCheckMKL(real matrix y, real scalar dir, real matrix X, real scalar d, real rowvector beta, real matrix Z,
	real rowvector opt, real matrix E, | real rowvector fm, real matrix tr)
{
	class LsmNat scalar N
	real colvector r
	real rowvector prior, pvar
	real matrix EPrior
	real scalar n, P

	n = rows(y)
	pvar = lsm_priorvar(X, n, dir, 9)
	prior = (1, 1, 1)
	EPrior = (opt[5], opt[6], opt[5] \ opt[7], opt[8], opt[7] \ opt[9], opt[10], opt[9])
	P = n * d + cols(beta) + n * (opt[1] + opt[2] + opt[3])
	N.open()
	if (args() < 9) fm = (0, 0, 0, 0)
	if (args() < 10) tr = J(0, 0, .)
	N.setup(y, dir, X, d, prior, pvar, beta, Z, opt[1..4], EPrior, E, J(1, 5, 0), fm, tr)
	r = N.call(sprintf("4 %f", N.L), 1 + 2 * P)
	N.close()
	return(r')
}

// plot data for nwlsmplot: node rows (x, y, colour group, weight, label) and edge rows (x1, y1, x2, y2) kept in externals until LsmPlotStore()
void LsmPlotPrep(string scalar Aname, string scalar Zname, real scalar d1, real scalar d2, real scalar dir, string scalar color, string scalar clsmat,
	string scalar sizevar, string scalar sizemat, string scalar label, real scalar edges)
{
	real matrix A, Z, ed, pd
	real colvector w, cg
	real scalar n, i, j
	string colvector lb
	external __nwlsm_pd, __nwlsm_ed, __nwlsm_lb

	A = st_matrix(Aname)
	Z = st_matrix(Zname)
	n = rows(Z)
	w = J(n, 1, 1)
	if (sizevar != "") w = st_data(., sizevar)
	if (sizemat != "") w = st_matrix(sizemat)[., 1]
	if (sizevar != "" | sizemat != "") {
		// marker area proportional to a positive weight: shift to a positive range
		w = w :- min(w)
		w = w :+ 0.15 * max((max(w), 1e-6)) :+ 1e-6
	}
	cg = J(n, 1, 1)
	if (color != "") cg = st_data(., color)
	if (clsmat != "") cg = st_matrix(clsmat)[., 1]
	lb = J(n, 1, "")
	if (label != "") lb = st_sdata(., label)
	pd = (Z[., d1], Z[., d2], cg, w)
	ed = J(0, 4, .)
	if (edges) {
		for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
			if (i == j) continue
			if (!dir & j < i) continue
			if (A[i, j] < . & A[i, j] != 0) ed = ed \ (Z[i, d1], Z[i, d2], Z[j, d1], Z[j, d2])
		}
	}
	__nwlsm_pd = pd
	__nwlsm_ed = ed
	__nwlsm_lb = lb
}

void LsmPlotStore()
{
	real scalar n, m
	external __nwlsm_pd, __nwlsm_ed, __nwlsm_lb

	n = rows(__nwlsm_pd)
	m = rows(__nwlsm_ed)
	st_store((1::n), ("x", "y", "cg", "w"), __nwlsm_pd)
	st_sstore((1::n), "lb", __nwlsm_lb)
	if (m > 0) st_store((n + 1::n + m), ("ex1", "ey1", "ex2", "ey2"), __nwlsm_ed)
}

// ---- BIC of the latent position model as latentnet's bic.ergmm (latentnet 2.12.0) ----
// log-likelihood of the Bernoulli tie model at the maximum over the coefficients for fixed positions and effects (Newton with step halving);
// y: n x n, X: stacked covariates (n nc x n), off: n x n offset (minus the distance plus the effects). Returns (lpY, beta).
real rowvector lsm_condmle(real matrix y, real scalar dir, real matrix X, real matrix off)
{
	real scalar n, nc, K, i, j, m, k, it, ll, ll1, st, ok
	real colvector yv, ov, eta, p, beta, g, nb, ho
	real matrix Xd, H

	n = rows(y)
	nc = rows(X) / n
	K = nc + 1
	yv = J(0, 1, .)
	ov = J(0, 1, .)
	Xd = J(0, K, .)
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j | (!dir & j < i) | y[i, j] >= .) continue
		yv = yv \ y[i, j]
		ov = ov \ off[i, j]
		Xd = Xd \ (1, (nc > 0 ? X[i :+ n * (0..(nc - 1)), j]' : J(1, 0, .)))
	}
	m = rows(yv)
	beta = J(K, 1, 0)
	eta = ov
	ll = sum(yv :* eta - (eta :* (eta :> 0) + ln(1 :+ exp(-abs(eta)))))
	for (it = 1; it <= 100; it++) {
		p = 1 :/ (1 :+ exp(-eta))
		g = cross(Xd, yv - p)
		H = cross(Xd, p :* (1 :- p), Xd)
		ho = cholinv(H + 1e-12 * I(K)) * g
		st = 1
		ok = 0
		for (k = 1; k <= 30; k++) {
			nb = beta + st * ho
			eta = ov + Xd * nb
			ll1 = sum(yv :* eta - (eta :* (eta :> 0) + ln(1 :+ exp(-abs(eta)))))
			if (ll1 >= ll - 1e-12) {
				ok = 1
				break
			}
			st = st / 2
		}
		if (!ok) break
		beta = nb
		if (abs(ll1 - ll) < 1e-12 * (1 + abs(ll))) {
			ll = ll1
			break
		}
		ll = ll1
	}
	return((ll, beta'))
}

// k-means (Lloyd, nstart random starts from data points); returns the cluster of every row, the centres, the within sums of squares
void lsm_kmeans(real matrix Z, real scalar G, real scalar nstart, real colvector cl, real matrix C, real colvector wss)
{
	real scalar n, d, s, it, g, best, i
	real colvector cl0, cl1, wss1, idx, w
	real matrix C1, DM

	n = rows(Z)
	d = cols(Z)
	best = .
	for (s = 1; s <= nstart; s++) {
		idx = order(runiform(n, 1), 1)[1::G]
		C1 = Z[idx, .]
		cl1 = J(n, 1, 0)
		for (it = 1; it <= 100; it++) {
			DM = J(n, G, 0)
			for (g = 1; g <= G; g++) DM[., g] = rowsum((Z :- C1[g, .]) :^ 2)
			cl0 = cl1
			cl1 = J(n, 1, .)
			for (i = 1; i <= n; i++) {
				minindex(DM[i, .], 1, w, idx)
				cl1[i] = w[1]
			}
			if (all(cl1 :== cl0)) break
			for (g = 1; g <= G; g++) if (any(cl1 :== g)) C1[g, .] = mean(select(Z, cl1 :== g))
		}
		wss1 = J(G, 1, 0)
		for (g = 1; g <= G; g++) wss1[g] = sum((select(Z, cl1 :== g) :- C1[g, .]) :^ 2)
		if (best >= . | sum(wss1) < best) {
			best = sum(wss1)
			cl = cl1
			C = C1
			wss = wss1
		}
	}
}

// EM for a mixture of G spherical normals with their own variances (latentnet's mbc.VII.EM); returns the log likelihood (missing if the fit failed)
real scalar lsm_em_vii(real matrix Z, real scalar G)
{
	real scalar n, d, att, it, g, dm, llk, tol, ok
	real colvector pK, v, cl, wss, sz, tot, pK1, v1
	real matrix mu, pz, lg, mu1
	transmorphic rs

	n = rows(Z)
	d = cols(Z)
	tol = sqrt(epsilon(1))
	rs = rngstate()
	rseed(20261119)
	llk = .
	for (att = 1; att <= 15; att++) {
		if (G > 1) {
			lsm_kmeans(Z, G, 15, cl, mu, wss)
			sz = J(G, 1, 0)
			for (g = 1; g <= G; g++) sz[g] = sum(cl :== g)
			pK = sz / n
			v = wss :/ (sz + (sz :== 0))
			for (g = 1; g <= G; g++) if (v[g] <= tol) v[g] = sum((Z :- mean(Z)) :^ 2) / (n - d)
		}
		else {
			mu = mean(Z)
			pK = 1
			v = sum((Z :- mean(Z)) :^ 2) / (n - d)
		}
		ok = 1
		for (it = 1; it <= 200; it++) {
			lg = J(n, G, .)
			for (g = 1; g <= G; g++) lg[., g] = ln(pK[g] + (pK[g] == 0) * 1e-300) :- 0.5 * d * ln(2 * pi() * v[g]) :- rowsum((Z :- mu[g, .]) :^ 2) / (2 * v[g])
			tot = rowmax(lg)
			pz = exp(lg :- tot)
			pz = pz :/ rowsum(pz)
			pK1 = mean(pz)'
			mu1 = J(G, d, 0)
			v1 = J(G, 1, 0)
			for (g = 1; g <= G; g++) {
				mu1[g, .] = colsum(pz[., g] :* Z) / sum(pz[., g])
				v1[g] = sum(pz[., g] :* rowsum((Z :- mu1[g, .]) :^ 2) / d) / sum(pz[., g])
			}
			if (hasmissing(mu1) | hasmissing(v1) | min(v1) <= 0) {
				ok = 0
				break
			}
			dm = max((mean(abs(pK1 - pK)) / mean(abs(pK)), mean(abs(vec(mu1 - mu))) / max((mean(abs(vec(mu))), 1e-300)), mean(abs(v1 - v)) / mean(abs(v))))
			pK = pK1
			mu = mu1
			v = v1
			if (dm < tol) break
		}
		if (!ok) continue
		if (max(v) / min(v) > 1 / sqrt(epsilon(1))) continue
		lg = J(n, G, .)
		for (g = 1; g <= G; g++) lg[., g] = ln(pK[g]) :- 0.5 * d * ln(2 * pi() * v[g]) :- rowsum((Z :- mu[g, .]) :^ 2) / (2 * v[g])
		tot = rowmax(lg)
		llk = sum(tot + ln(rowsum(exp(lg :- tot))))
		break
	}
	rngstate(rs)
	return(llk)
}

// BIC components of latentnet for the positions Z (n x d), effects E (n x 3: sender, receiver, sociality), G clusters (0: none), eff the effective
// number of observations; returns (Y, Z, sender, receiver, sociality, overall, p, lpY, G used)
real rowvector LsmBic(real matrix y, real scalar dir, real matrix X, real matrix Z, real matrix E, real rowvector act, real scalar G, real scalar eff)
{
	real scalar n, d, i, j, p, lpY, Gs, llk, s2, bY, bZ, t
	real matrix off
	real rowvector out, r

	n = rows(y)
	d = cols(Z)
	off = J(n, n, 0)
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j) continue
		off[i, j] = -sqrt(sum((Z[i, .] - Z[j, .]) :^ 2))
		if (act[1]) off[i, j] = off[i, j] + E[i, 1]
		if (act[2]) off[i, j] = off[i, j] + E[j, 2]
		if (act[3]) off[i, j] = off[i, j] + E[i, 3] + E[j, 3]
	}
	r = lsm_condmle(y, dir, X, off)
	lpY = r[1]
	p = cols(r) - 1
	bY = -2 * lpY + p * ln(eff)
	Gs = G
	if (G > 0) {
		llk = .
		while (llk >= . & Gs > 0) {
			llk = lsm_em_vii(Z, Gs)
			if (llk >= .) Gs--
		}
		if (Gs < G) errprintf("nwlsm: bad clustering, treating %g clusters as empty\n", G - Gs)
		bZ = -2 * llk + (Gs - 1 + d * Gs + Gs) * ln(n)
	}
	else {
		s2 = sum(Z :^ 2) / n
		bZ = n * d * ln(2 * pi() * s2) + n + ln(n * d)
	}
	out = (bY, bZ, 0, 0, 0)
	for (t = 1; t <= 3; t++) if (act[t]) {
		s2 = sum(E[., t] :^ 2) / n
		out[2 + t] = n * (ln(2 * pi() * s2) + 1) + ln(n)
	}
	return((out, sum(out), p, lpY, Gs))
}

// BIC of the last fit (positions: the MKL positions; effects E: n x 3 as e(mkl_effects)); eff: effective number of observations (0: the number of ties)
real rowvector LsmBicFit(real matrix E, real scalar G, real scalar eff)
{
	external __nwlsm_gY, __nwlsm_gX, __nwlsm_gm, __nwlsm_MZ
	real matrix y, A
	real scalar dir, n

	y = __nwlsm_gY
	n = rows(y)
	dir = __nwlsm_gm[1]
	if (eff <= 0) {
		A = (y :== 1) :* (J(n, n, 1) - I(n))
		eff = (dir ? sum(A) : sum(lowertriangle(A, 0)))
	}
	return((LsmBic(y, dir, __nwlsm_gX, __nwlsm_MZ, E, __nwlsm_gm[3..5], G, eff), eff))
}

end
