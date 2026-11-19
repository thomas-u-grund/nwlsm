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
capture mata: mata drop lsm_ztrans()
capture mata: mata drop LsmPlotPrep()
capture mata: mata drop LsmPlotStore()

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

real scalar LSM_PROTOCOL_NEEDED() return(3)

class LsmNat {
	real scalar isopen, L
	string scalar origframe
	void open()
	void close()
	void setup()
	real colvector call()
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

// y: n x n (0/1, missing = .), X: (nc*n) x n, the covariate matrices stacked, d: dimension, prior: (scale, df, s2 start),
// pvar, b0: 1 x K, Z0: n x d, flags: (sender, receiver, sociality, bilinear), eprior: 3 x 3 (scale, df, start variance of the effects),
// E0: n x 3 start values of the effects
void LsmNat::setup(real matrix y, real scalar directed, real matrix X, real scalar d, real rowvector prior,
		real rowvector pvar, real rowvector b0, real matrix Z0, real rowvector flags, real matrix eprior, real matrix E0, real rowvector cl)
{
	real scalar n, nc, k
	real colvector blob

	n = rows(y)
	nc = rows(X) / n
	blob = (LSM_PROTOCOL_NEEDED() \ n \ directed \ nc \ d \ (flags[1] + 2 * flags[2] + 4 * flags[3]) \ flags[4] \ cl[1] \ cl[2] \ cl[3] \ cl[4] \ cl[5] \ J(4, 1, 0)) \ prior' \ pvar' \ b0' \ vec(Z0') \ vec(eprior') \ vec(E0) \ vec(editmissing(y, -1)')
	if (nc > 0) {
		for (k = 1; k <= nc; k++) blob = blob \ vec(editmissing(X[((k - 1) * n + 1)..(k * n), .], 0)')
	}
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
	real rowvector opt, string rowvector nm)
{
	class LsmNat scalar N
	real scalar n, nc, K, S, Sc, C, G, ng, k, it, rows_, off, nd, matok, na, t, i, g2, bi
	real rowvector prior, pvar, b0, flags, mb, act, cl
	real matrix Z0, beta, Zs, Dbar, Zm, Za, Zr, T, EP, Pm, VT, Es, Ev, Eall, EPrior, E0, Em, S2, RH, Mu, Lam, Kd, PI, Q, S2r, Lr, Kr, MUa, R, Zi
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
	lsm_start(y, dir, X, d, random, seed, Z0, b0, flags[4])
	N.open()
	N.setup(y, dir, X, d, prior, pvar, b0, Z0, flags, EPrior, E0, cl)
	rows_ = 2 * K + 8 + S * K + S + S * ng + 3 * S + S * nd + na * S * n + (G > 0 ? S * G * d + S * G + S * n : 0) + 2 * n * n + nd + K + na * n + 3
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
	T = J(K, 5, .)
	T[., 1] = mean(beta)'
	T[., 2] = sqrt(diagonal(variance(beta)))
	T[., 3] = lsm_mcse(beta, C)'
	for (k = 1; k <= K; k++) {
		srt = sort(beta[., k], 1)
		T[k, 4] = srt[max((1, ceil(0.025 * S)))]
		T[k, 5] = srt[min((S, floor(0.975 * S) + 1))]
	}
	st_matrix(nm[1], mean(beta))
	st_matrix(nm[2], variance(beta))
	st_matrix(nm[3], T)
	// variance parameters: sigma2 of the positions (one per cluster with clusters), then the active effects
	VT = J(0, 5, .)
	for (k = 1; k <= ng; k++) VT = VT \ lsm_vartable(S2[., k], C)
	for (t = 1; t <= 3; t++) if (act[t]) VT = VT \ lsm_vartable(Ev[., t], C)
	st_matrix(nm[11], VT)
	// Gelman-Rubin statistics (several chains): coefficients, log likelihood, variance parameters
	if (C > 1) {
		RH = lsm_rhat(beta, C)' \ lsm_rhat(llk, C)
		for (k = 1; k <= ng; k++) RH = RH \ lsm_rhat(S2[., k], C)
		for (t = 1; t <= 3; t++) if (act[t]) RH = RH \ lsm_rhat(Ev[., t], C)
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
	st_local("__priorscale", strofreal(G > 0 ? cl[2] : prior[1], "%18.10g"))
	st_local("__priordf", strofreal(G > 0 ? cl[3] : prior[2], "%18.10g"))
	st_matrix(nm[4] + "_pv", pvar)
	st_matrix(nm[4] + "_b0", b0)
	external __nwlsm_beta, __nwlsm_llk, __nwlsm_s2, __nwlsm_Z, __nwlsm_D, __nwlsm_ZA, __nwlsm_ev, __nwlsm_E, __nwlsm_P, __nwlsm_MZ, __nwlsm_C, __nwlsm_K, __nwlsm_lam
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
	external __nwlsm_beta, __nwlsm_llk, __nwlsm_s2, __nwlsm_Z, __nwlsm_D, __nwlsm_ZA, __nwlsm_ev, __nwlsm_E, __nwlsm_P, __nwlsm_MZ, __nwlsm_C, __nwlsm_K, __nwlsm_lam

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
void LsmSaveDraws(string rowvector names, real rowvector act)
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
}

// test helper: the plugin's log likelihood and log prior of Z at (beta, Z, s2), and its incremental log-likelihood change
// if actor i moves to znew. Returns (loglik, logprior_Z, dll, dll_effect), the last for the effect of type t (0 = none) of actor i
// changing by delta. E: n x 3 start values of the effects; opt as for LsmFit.
real rowvector LsmCheck(real matrix y, real scalar dir, real matrix X, real scalar d, real rowvector beta, real matrix Z,
	real scalar s2, real scalar i, real rowvector znew, real rowvector opt, real matrix E, real scalar t, real scalar delta)
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
	N.setup(y, dir, X, d, prior, pvar, beta, Z, opt[1..4], EPrior, E, J(1, 5, 0))
	r1 = N.call(sprintf("1 %f", N.L), 2)
	r2 = N.call(sprintf("2 %f %f ", N.L, i) + invtokens(strofreal(znew, "%20.15g")), 1)
	r3 = 0
	if (t > 0) r3 = N.call(sprintf("3 %f %f %f %s", N.L, t - 1, i, strofreal(delta, "%20.15g")), 1)
	N.close()
	return((r1', r2, r3))
}

// test helper: objective of the MKL optimisation at the state (target = 0.9 for a tie, 0.15 for none), the analytic gradient and a
// central-difference gradient. Returns (f, analytic' , numeric') as a row vector.
real rowvector LsmCheckMKL(real matrix y, real scalar dir, real matrix X, real scalar d, real rowvector beta, real matrix Z,
	real rowvector opt, real matrix E)
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
	N.setup(y, dir, X, d, prior, pvar, beta, Z, opt[1..4], EPrior, E, J(1, 5, 0))
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

end
