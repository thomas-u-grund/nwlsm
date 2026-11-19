* Shared by the influence-model tests: lsmi_sim simulates data from the exact model of Sweet and Adhikari (2020) with known parameters
* and lsmi_rep fits it and records the result. Included by do-files after nwlsm_testlib.do.
capture program drop lsmi_sim
* network: logit P(A_ij = 1) = a0 - |Z_i - Z_j| (Z_i ~ N(0, zsd^2 I_d)); outcome: y = b0 + b1 ylag + b3 x + b2 m + e,  e ~ N(0, s2),
* m_k = weighted mean of the lagged outcomes of the G nearest actors (G = 0: all), weights 1/d (wt 0) or exp(-d) (wt 1). Variables: __a1..__an
* (network), y, ylag, x in the n observations; the true distances, positions and weights are kept in Mata externals.
program lsmi_sim
	syntax, N(integer) DIR(integer) D(integer) A0(real) B0(real) B1(real) B2(real) B3(real) S2(real) G(integer) WT(integer) SEED(integer) ZSD(real)
	clear
	set obs `n'
	forvalues j = 1/`n' {
		quietly gen double __a`j' = .
	}
	quietly gen double y = .
	quietly gen double ylag = .
	quietly gen double x = .
	mata: lsmi_simulate(`n', `dir', `d', `a0', `b0', `b1', `b2', `b3', `s2', `g', `wt', `seed', `zsd')
end
mata:
mata set matastrict off
void lsmi_simulate(real scalar n, real scalar dir, real scalar d, real scalar a0, real scalar b0, real scalar b1, real scalar b2, real scalar b3,
	real scalar s2, real scalar G, real scalar wt, real scalar seed, real scalar zsd)
{
	real matrix Z, D, A, W
	real colvector ylag, x, y, o, w
	real scalar k, Gm
	external real matrix __lsmi_Z, __lsmi_D, __lsmi_W

	rseed(seed)
	Z = rnormal(n, d, 0, zsd)
	D = lsmt_pd(Z)
	A = (runiform(n, n) :< invlogit(a0 :- D))
	if (!dir) {
		A = uppertriangle(A, 0)
		A = A + A'
	}
	_diag(A, 0)
	ylag = rnormal(n, 1, 0, 1)
	x = rnormal(n, 1, 0, 1)
	Gm = (G <= 0 | G >= n - 1) ? n - 1 : G
	W = J(n, n, 0)
	for (k = 1; k <= n; k++) {
		o = order(D[k, .]', 1)
		o = select(o, o :!= k)
		o = o[1..Gm]
		w = wt == 0 ? 1 :/ D[k, o]' : exp(-D[k, o]')
		W[k, o] = (w / sum(w))'
	}
	y = b0 :+ b1 * ylag + b3 * x + b2 * (W * ylag) + rnormal(n, 1, 0, sqrt(s2))
	st_store(., (1..n), A)
	st_store(., n + 1, y)
	st_store(., n + 2, ylag)
	st_store(., n + 3, x)
	__lsmi_Z = Z; __lsmi_D = D; __lsmi_W = W
}
// correlation of the upper triangles of two square matrices
real scalar lsmi_cor_ut(real matrix A, real matrix B)
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
// Procrustes root mean square distance (positions centred, optimal rotation or reflection, no scaling)
real scalar lsmi_procrmse(real matrix Z, real matrix R)
{
	real matrix A, B, U, V, Rot
	real rowvector s
	A = Z :- mean(Z)
	B = R :- mean(R)
	svd(A' * B, U, s, V)
	Rot = U * V'
	return(sqrt(sum((A * Rot - B) :^ 2) / rows(A)))
}
end
