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
capture mata: mata drop lsm_procrustes()
capture mata: mata drop lsm_ymat()
capture mata: mata drop lsm_node_covs()
capture mata: mata drop LsmFit()
capture mata: mata drop LsmCheck()
capture mata: mata drop LsmSample()
capture mata: mata drop LsmSaveDraws()

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

real scalar LSM_PROTOCOL_NEEDED() return(1)

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
// pvar, b0: 1 x K, Z0: n x d
void LsmNat::setup(real matrix y, real scalar directed, real matrix X, real scalar d, real rowvector prior,
		real rowvector pvar, real rowvector b0, real matrix Z0)
{
	real scalar n, nc
	real colvector blob

	n = rows(y)
	nc = rows(X) / n
	blob = (LSM_PROTOCOL_NEEDED() \ n \ directed \ nc \ d \ J(11, 1, 0)) \ prior' \ pvar' \ b0' \ vec(Z0') \ vec(editmissing(y, -1)') 
	if (nc > 0) {
		real scalar k
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
	real matrix Z0, real rowvector b0)
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
	off = -vec(sqrt(Dm))[r]
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

// Monte Carlo standard error of the mean of each column (batch means, batches of about sqrt(S) draws)
real rowvector lsm_mcse(real matrix x)
{
	real scalar S, b, nb, k
	real matrix bm

	S = rows(x)
	b = max((1, floor(sqrt(S))))
	nb = floor(S / b)
	if (nb < 2) return(sqrt(diagonal(variance(x))' / S))
	bm = J(nb, cols(x), .)
	for (k = 1; k <= nb; k++) bm[k, .] = mean(x[((k - 1) * b + 1)..(k * b), .])
	return(sqrt(diagonal(variance(bm))' / nb))
}

// Procrustes (rotation and reflection, no scaling) of Z onto Ref, both centred; returns the aligned Z
real matrix lsm_procrustes(real matrix Z, real matrix Ref)
{
	real matrix M, U, Vt
	real rowvector s

	M = cross(Z, Ref)
	svd(M, U, s, Vt)
	return(Z * (U * Vt))
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

/*
	LsmFit: runs the chain and stores the summary. Posterior sample stays in the Mata externals __nwlsm_beta (S x K), __nwlsm_llk,
	__nwlsm_s2 (S x 1), __nwlsm_Z (S x n d, by actor then dimension) and __nwlsm_D (the mean distance matrix).
	Results into Stata matrices: B (1 x K means), V (K x K covariance), T (K x 5: mean, sd, mcse, 2.5%, 97.5%), D (n x n, if it fits),
	Z (n x d aligned mean positions, if it fits), ACC (1 x (K+1) acceptance rates, step sizes below it).
*/
void LsmFit(real matrix y, real scalar dir, real matrix X, real scalar d, real scalar nsample, real scalar burnin,
	real scalar interval, real scalar seed, real scalar pilot, real scalar plen, real scalar random,
	string scalar Bn, string scalar Vn, string scalar Tn, string scalar Dn, string scalar Zn, string scalar An)
{
	class LsmNat scalar N
	real scalar n, nc, K, S, k, it, rows_, off, refdraw, nd, matok, lp, bestlp
	real rowvector prior, pvar, b0, q, lo, hi
	real matrix Z0, beta, Zs, Dbar, Zm, Za, Zr, T, ZZ
	real colvector res, llk, s2, srt

	n = rows(y)
	nc = rows(X) / n
	K = nc + 1
	S = nsample
	nd = n * d
	pvar = lsm_priorvar(X, n, dir, 9)
	prior = (0.125 * n^(2 / d), sqrt(n), 0.125 * n^(2 / d))
	lsm_start(y, dir, X, d, random, seed, Z0, b0)
	N.open()
	N.setup(y, dir, X, d, prior, pvar, b0, Z0)
	rows_ = 2 * K + 2 + S * K + 2 * S + S * nd + n * n
	res = N.call(sprintf("10 %f %f %f %f %f %f %f %f %f %f", N.L, seed, burnin, interval, S, pilot, plen, 0.6, 0.6, 0.234), rows_)
	N.close()
	if (hasmissing(res)) {
		errprintf("nwlsm: the sampler returned missing values\n")
		exit(498)
	}
	beta = colshape(res[(2 * K + 3)..(2 * K + 2 + S * K)], K)
	off = 2 * K + 2 + S * K
	llk = res[(off + 1)..(off + S)]
	s2 = res[(off + S + 1)..(off + 2 * S)]
	off = off + 2 * S
	Zs = colshape(res[(off + 1)..(off + S * nd)], nd)
	off = off + S * nd
	Dbar = colshape(res[(off + 1)..(off + n * n)], n)
	// summary of the coefficients
	T = J(K, 5, .)
	T[., 1] = mean(beta)'
	T[., 2] = sqrt(diagonal(variance(beta)))
	T[., 3] = lsm_mcse(beta)'
	for (k = 1; k <= K; k++) {
		srt = sort(beta[., k], 1)
		T[k, 4] = srt[max((1, ceil(0.025 * S)))]
		T[k, 5] = srt[min((S, floor(0.975 * S) + 1))]
	}
	st_matrix(Bn, mean(beta))
	st_matrix(Vn, variance(beta))
	st_matrix(Tn, T)
	// positions: centre, align every draw to the reference (the draw with the highest log posterior), then refine twice
	bestlp = .
	refdraw = 1
	for (it = 1; it <= S; it++) {
		lp = llk[it] - sum(Zs[it, .] :^ 2) / (2 * s2[it]) - 0.5 * nd * ln(s2[it])
		if (bestlp == . | lp > bestlp) {
			bestlp = lp
			refdraw = it
		}
	}
	Zr = colshape(Zs[refdraw, .], d)
	Zr = Zr :- mean(Zr)
	for (k = 1; k <= 3; k++) {
		Zm = J(n, d, 0)
		for (it = 1; it <= S; it++) {
			Za = colshape(Zs[it, .], d)
			Za = Za :- mean(Za)
			Zm = Zm + lsm_procrustes(Za, Zr)
		}
		Zm = Zm / S
		Zr = Zm :- mean(Zm)
	}
	matok = (n <= c("max_matdim"))
	if (matok) {
		st_matrix(Dn, Dbar)
		st_matrix(Zn, Zr)
	}
	st_matrix(An, (res[1..(K + 1)]' \ res[(K + 2)..(2 * K + 2)]'))
	st_local("__matok", strofreal(matok))
	st_local("__llk", strofreal(mean(llk), "%18.10g"))
	st_local("__llksd", strofreal(sqrt(variance(llk)), "%18.10g"))
	st_local("__s2", strofreal(mean(s2), "%18.10g"))
	st_local("__s2sd", strofreal(sqrt(variance(s2)), "%18.10g"))
	st_local("__priorscale", strofreal(prior[1], "%18.10g"))
	st_local("__priordf", strofreal(prior[2], "%18.10g"))
	st_matrix(Dn + "_pv", pvar)
	st_matrix(Dn + "_b0", b0)
	external __nwlsm_beta, __nwlsm_llk, __nwlsm_s2, __nwlsm_Z, __nwlsm_D, __nwlsm_ZA
	__nwlsm_beta = beta
	__nwlsm_llk = llk
	__nwlsm_s2 = s2
	__nwlsm_Z = Zs
	__nwlsm_D = Dbar
	__nwlsm_ZA = Zr
}

// the posterior sample kept in Mata by the last fit: LsmSample("beta"), "llk", "s2", "Z" (one row per draw), "D", "Zalign"
transmorphic matrix LsmSample(string scalar what)
{
	external __nwlsm_beta, __nwlsm_llk, __nwlsm_s2, __nwlsm_Z, __nwlsm_D, __nwlsm_ZA

	if (what == "beta") return(__nwlsm_beta)
	if (what == "llk") return(__nwlsm_llk)
	if (what == "s2") return(__nwlsm_s2)
	if (what == "Z") return(__nwlsm_Z)
	if (what == "D") return(__nwlsm_D)
	if (what == "Zalign") return(__nwlsm_ZA)
	errprintf("LsmSample(): unknown item\n")
	exit(198)
}

// writes the draws of the coefficients, the log likelihood and the latent-position variance into the (empty) current dataset
void LsmSaveDraws(string rowvector names)
{
	real scalar k, S
	external __nwlsm_beta, __nwlsm_llk, __nwlsm_s2

	S = rows(__nwlsm_beta)
	st_addobs(S)
	for (k = 1; k <= cols(__nwlsm_beta); k++) {
		(void) st_addvar("double", names[k])
		st_store(., names[k], __nwlsm_beta[., k])
	}
	(void) st_addvar("double", "loglik")
	st_store(., "loglik", __nwlsm_llk)
	(void) st_addvar("double", "sigma2")
	st_store(., "sigma2", __nwlsm_s2)
}

// test helper: the plugin's log likelihood and log prior of Z at (beta, Z, s2), and its incremental log-likelihood change
// if actor i moves to znew. Returns (loglik, logprior_Z, dll).
real rowvector LsmCheck(real matrix y, real scalar dir, real matrix X, real scalar d, real rowvector beta, real matrix Z,
	real scalar s2, real scalar i, real rowvector znew)
{
	class LsmNat scalar N
	real colvector r1, r2
	real rowvector prior, pvar
	real scalar n

	n = rows(y)
	pvar = lsm_priorvar(X, n, dir, 9)
	prior = (1, 1, s2)
	N.open()
	N.setup(y, dir, X, d, prior, pvar, beta, Z)
	r1 = N.call(sprintf("1 %f", N.L), 2)
	r2 = N.call(sprintf("2 %f %f ", N.L, i) + invtokens(strofreal(znew, "%20.15g")), 1)
	N.close()
	return((r1', r2))
}

end
