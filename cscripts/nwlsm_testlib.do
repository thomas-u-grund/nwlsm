* Shared set-up of the nwlsm tests: scratch PLUS directory, ad-paths, a network from an edge list with actors 1..n.
* Run the tests from the root of the repository. NWCORE: directory of the core package (default: the owner's checkout).
capture program drop nwinstall
local core : env NWCORE
if "`core'" == "" local core "/Users/tgrund/FILES_NEW/SOFTWARE/nwcommands"
local scratch : env NWLSM_SCRATCH
if "`scratch'" == "" local scratch "/private/tmp/claude-501/-Users-tgrund-FILES-NEW-SOFTWARE-nwcommands/9d5243e9-ee49-4558-835d-02925b410730/scratchpad/nwlsm_work"
local tag : env NWLSM_TAG
tempfile __f
local plus "`scratch'/plus_`tag'_`=floor(runiform()*1e9)'"
capture mkdir "`scratch'"
capture mkdir "`plus'"
sysdir set PLUS "`plus'"
local repo "`c(pwd)'"
adopath ++ "`core'"
adopath ++ "`core'/lib"
adopath ++ "`repo'"
adopath ++ "`repo'/lib"
adopath ++ "`repo'/lib/plugins/macos"
global LSM_REPO "`repo'"
global LSM_SCRATCH "`scratch'"

capture program drop lsmt_mknet
* lsmt_mknet name, edges(csv) n(#) dir(0|1): network of n actors, actor k = row k, from an edge list (from, to) with ids 1..n
program lsmt_mknet
	syntax name(name=nm), EDGES(string) N(integer) DIR(integer)
	preserve
	import delimited `"`edges'"', clear varnames(1)
	mata: lsmt_edges(`n', `dir')
	restore
	preserve
	clear
	set obs `n'
	forvalues j = 1/`n' {
		quietly gen double __a`j' = .
	}
	mata: st_store(., ., st_matrix("__lsmt_A"))
	if `dir' nwset __a*, name(`nm') directed
	else nwset __a*, name(`nm') undirected
	restore
	capture matrix drop __lsmt_A
end

mata:
mata set matastrict off
void lsmt_edges(real scalar n, real scalar dir)
{
	real matrix A, E
	real scalar k

	E = st_data(., ("from", "to"))
	A = J(n, n, 0)
	for (k = 1; k <= rows(E); k++) {
		A[E[k, 1], E[k, 2]] = 1
		if (!dir) A[E[k, 2], E[k, 1]] = 1
	}
	st_matrix("__lsmt_A", A)
}

// matrix of the Euclidean distances between the rows of Z
real matrix lsmt_pd(real matrix Z)
{
	real matrix D
	real scalar n, i, j

	n = rows(Z)
	D = J(n, n, 0)
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) D[i, j] = sqrt(sum((Z[i, .] - Z[j, .]) :^ 2))
	return(D)
}
// full log likelihood of the model by direct computation. E: n x 3 effects (sender, receiver, sociality), opt: (sender, receiver, sociality, bilinear, ...)
real scalar lsmt_ll(real matrix y, real scalar dir, real matrix X, real rowvector beta, real matrix Z, real matrix E, real rowvector opt)
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
		t = t + y[i, j] * eta - ln(1 + exp(eta))
	}
	return(t)
}
end
