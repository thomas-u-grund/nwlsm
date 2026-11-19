* Valued networks, more: (1) Poisson with two latent clusters against the planted clusters and latentnet; (2) results independent of threads() for the three
* families (chains, clusters, sender/receiver effects); (3) error paths of family(), trials(), dispvar() and valued networks; (4) the valued-network
* flag of the core: a valued network without family() is refused unless it holds a binary graph; family(poisson) on a binary network.
do cscripts/nwlsm_testlib.do
set linesize 255
mata:
mata set matastrict off
real scalar lsmt_ari(real colvector a, real colvector b)
{
	real matrix T
	real scalar n, i, sij, sa, sb, ex, mx
	real colvector ua, ub

	ua = uniqrows(a)
	ub = uniqrows(b)
	T = J(rows(ua), rows(ub), 0)
	for (i = 1; i <= rows(a); i++) T[selectindex(ua :== a[i]), selectindex(ub :== b[i])] = T[selectindex(ua :== a[i]), selectindex(ub :== b[i])] + 1
	n = rows(a)
	sij = sum(T :* (T :- 1) / 2)
	sa = sum(rowsum(T) :* (rowsum(T) :- 1) / 2)
	sb = sum(colsum(T)' :* (colsum(T)' :- 1) / 2)
	ex = sa * sb / (n * (n - 1) / 2)
	mx = (sa + sb) / 2
	if (mx == ex) return(1)
	return((sij - ex) / (mx - ex))
}
end
* (1) clusters
lsmt_mkvnet net_poisclu36, y("$LSM_REPO/dev/ref/poisclu36_y.csv") dir(0)
preserve
import delimited "$LSM_REPO/dev/ref/poisclu36_truth.csv", clear varnames(1)
mkmat truth, matrix(TR)
import delimited "$LSM_REPO/dev/ref/poisclu36_class.csv", clear varnames(1)
mkmat class, matrix(LC)
import delimited "$LSM_REPO/dev/ref/poisclu36_beta.csv", clear varnames(1)
mkmat mean sd, matrix(RB)
import delimited "$LSM_REPO/dev/ref/poisclu36_c2_beta.csv", clear varnames(1)
mkmat mean, matrix(RB2)
restore
local ok 1
matrix SB = J(4, 1, .)
matrix SM = J(4, 1, .)
forvalues s = 1/4 {
	quietly nwlsm net_poisclu36, family(poisson) dim(2) clusters(2) seed(`=3000 + `s'')
	matrix SB[`s', 1] = e(b)[1, 1]
	matrix SM[`s', 1] = e(summary)[1, 3]
	mata: st_local("a1", strofreal(lsmt_ari(st_matrix("e(class)"), st_matrix("TR")), "%7.4f"))
	mata: st_local("a2", strofreal(lsmt_ari(st_matrix("e(class)"), st_matrix("LC")), "%7.4f"))
	di as text "poisson clusters seed `s': adjusted Rand index with the planted clusters " `a1' ", with latentnet's classification " `a2'
	if `a1' < 0.9 | `a2' < 0.9 local ok 0
}
mata: st_local("avg", strofreal(mean(st_matrix("SB"))))
mata: st_local("se", strofreal(sqrt(mean(st_matrix("SM") :^ 2) / 4)))
local ref = (RB[1, 1] + RB2[1, 1]) / 2
local seref = max(abs(RB[1, 1] - RB2[1, 1]) / 2, RB[1, 2] / 7)
local z = abs(`avg' - `ref') / sqrt(`se'^2 + `seref'^2)
di as text "poisson clusters intercept: nwlsm " %8.4f `avg' "  latentnet " %8.4f `ref' "  diff/se = " %5.2f `z'
if `z' > 4 local ok 0
if !`ok' {
	di as error "FAIL poisson clusters"
	exit 9
}
di as result "PASS poisson clusters"
* (2) threads
capture program drop lsmt_samev
program lsmt_samev
	syntax anything(name=case), OPTS(string) T1(integer) T2(integer) NET(string)
	quietly nwlsm `net', `opts' threads(`t1')
	matrix b1 = e(b)
	matrix Z1 = e(Z)
	matrix P1 = e(P)
	matrix V1 = e(vsummary)
	local l1 = e(ll_mean)
	local m1 = e(mkl_obj)
	quietly nwlsm `net', `opts' threads(`t2')
	mata: st_local("db", strofreal(max(abs(st_matrix("b1") - st_matrix("e(b)")))))
	mata: st_local("dz", strofreal(max(abs(st_matrix("Z1") - st_matrix("e(Z)")))))
	mata: st_local("dp", strofreal(max(abs(st_matrix("P1") - st_matrix("e(P)")))))
	mata: st_local("dv", strofreal(max(abs(st_matrix("V1") - st_matrix("e(vsummary)")))))
	di as text "`case': threads(`t1') against threads(`t2'): max |diff| b " `db' "  Z " `dz' "  P " `dp' "  variances " `dv' "  ll " `=abs(`l1' - e(ll_mean))' "  mkl " `=abs(`m1' - e(mkl_obj))'
	if `db' != 0 | `dz' != 0 | `dp' != 0 | `dv' != 0 | `l1' != e(ll_mean) | `m1' != e(mkl_obj) {
		di as error "FAIL `case': the result depends on the number of threads"
		exit 9
	}
	di as result "PASS `case'"
end
lsmt_mkvnet net_pois30, y("$LSM_REPO/dev/ref/pois30_y.csv") dir(1)
lsmt_mkvnet net_bin30, y("$LSM_REPO/dev/ref/bin30_y.csv") dir(1)
lsmt_mkvnet net_gau30, y("$LSM_REPO/dev/ref/gau30_y.csv") dir(1)
lsmt_samev poisson, net(net_pois30) opts(family(poisson) dim(2) seed(5) nsample(300) burnin(500) chains(3) sender receiver) t1(1) t2(3)
lsmt_samev binomial, net(net_bin30) opts(family(binomial) trials(6) dim(2) seed(6) nsample(300) burnin(500) chains(2) sociality mklstarts(3)) t1(1) t2(4)
lsmt_samev normal, net(net_gau30) opts(family(normal) dispvar(1) dispdf(3) dim(2) seed(7) nsample(300) burnin(500) chains(3) bilinear) t1(1) t2(2)
lsmt_samev poisclu, net(net_poisclu36) opts(family(poisson) dim(2) seed(8) nsample(300) burnin(500) chains(2) clusters(2)) t1(1) t2(2)
* (3) errors
capture nwlsm net_pois30, family(gamma)
assert _rc == 198
capture nwlsm net_pois30, family(binomial)
assert _rc == 198
capture nwlsm net_pois30, family(poisson) trials(5)
assert _rc == 198
capture nwlsm net_pois30, family(normal) dispvar(0)
assert _rc == 198
capture nwlsm net_pois30, family(normal) dispdf(-1)
assert _rc == 198
capture nwlsm net_pois30, family(binomial) trials(0)
assert _rc == 198
capture nwlsm net_pois30, family(binomial) trials(2.5)
assert _rc == 198
capture nwlsm net_pois30, family(binomial) trials(2) nsample(50) burnin(10)
assert _rc == 198
capture nwlsm net_pois30, family(binomial) trials(nosuchthing)
assert _rc == 198
mata: st_matrix("__BADT", J(5, 5, 3))
capture nwlsm net_pois30, family(binomial) trials(__BADT)
assert _rc == 198
* a valued network needs family() unless it holds a binary graph
capture nwlsm net_pois30
assert _rc == 198
capture nwlsm net_gau30, nsample(50) burnin(10)
assert _rc == 198
* Poisson needs counts: the normal data are refused
capture nwlsm net_gau30, family(poisson)
assert _rc == 198
* an undirected valued network: sender is refused as for binary networks
lsmt_mkvnet net_pois30s, y("$LSM_REPO/dev/ref/pois30s_y.csv") dir(0)
capture nwlsm net_pois30s, family(poisson) sender
assert _rc == 198
* (4) family(bernoulli) is the default; a binary network with family(poisson) is fitted; the e() results
lsmt_mknet net40, edges("$LSM_REPO/dev/ref/random40_d2_edges.csv") n(40) dir(0)
quietly nwlsm net40, family(bernoulli) nsample(100) burnin(100) seed(3)
assert "`e(family)'" == "bernoulli"
quietly nwlsm net40, family(poisson) nsample(100) burnin(100) seed(3)
assert "`e(family)'" == "poisson"
quietly nwlsm net_gau30, family(normal) dispvar(2) dispdf(5) nsample(100) burnin(100) seed(3)
assert "`e(family)'" == "normal"
assert e(dispvar) == 2 & e(dispdf) == 5
matrix VT = e(vsummary)
local __rn : rownames VT
assert word("`__rn'", rowsof(VT)) == "dispersion"
quietly nwlsm net_bin30, family(binomial) trials(6) nsample(100) burnin(100) seed(3)
assert e(trials) == 6
* the draws of the dispersion are kept and saved
quietly nwlsm net_gau30, family(normal) dispvar(2) dispdf(5) nsample(100) burnin(100) seed(3) savedraws("$LSM_SCRATCH/drawsv.dta")
preserve
use "$LSM_SCRATCH/drawsv.dta", clear
confirm variable dispersion
assert _N == 100
summarize dispersion, meanonly
assert r(min) > 0
restore
display "PASS test_nwlsm_valued_misc: Poisson clusters, threads, error paths and e() results of the valued families"
