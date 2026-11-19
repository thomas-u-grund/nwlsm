* Threads: for a given seed and chains() the results must be identical for every threads() value (exactly, not up to Monte Carlo error).
*  (1) several chains: threads(1), 2 and 5 (more threads than chains);
*  (2) one chain on a network large enough for the inner per-dyad passes (n = 200: up to 4 threads work inside the chain and the MKL fit): threads(1) and threads(4);
*  (3) clusters and chains: threads(1) and threads(3);
*  (4) errors of the options.
do cscripts/nwlsm_testlib.do
capture program drop lsmt_same
program lsmt_same
	syntax anything(name=case), OPTS(string) T1(integer) T2(integer) NET(string)
	quietly nwlsm `net', `opts' threads(`t1')
	matrix b1 = e(b)
	matrix Z1 = e(Z)
	matrix D1 = e(D)
	matrix P1 = e(P)
	local l1 = e(ll_mean)
	local m1 = e(mkl_obj)
	quietly nwlsm `net', `opts' threads(`t2')
	mata: st_local("db", strofreal(max(abs(st_matrix("b1") - st_matrix("e(b)")))))
	mata: st_local("dz", strofreal(max(abs(st_matrix("Z1") - st_matrix("e(Z)")))))
	mata: st_local("dd", strofreal(max(abs(st_matrix("D1") - st_matrix("e(D)")))))
	mata: st_local("dp", strofreal(max(abs(st_matrix("P1") - st_matrix("e(P)")))))
	di as text "`case': threads(`t1') against threads(`t2'): max |diff| b " `db' "  Z " `dz' "  D " `dd' "  P " `dp' "  ll " `=abs(`l1' - e(ll_mean))' "  mkl " `=abs(`m1' - e(mkl_obj))'
	if `db' != 0 | `dz' != 0 | `dd' != 0 | `dp' != 0 | `l1' != e(ll_mean) | `m1' != e(mkl_obj) {
		di as error "FAIL `case': the result depends on the number of threads"
		exit 9
	}
	di as result "PASS `case'"
end
mata:
void lsmt_tie(real matrix A, real matrix Zt, real scalar i, real scalar j)
{
	if (runiform(1, 1) < invlogit(1 - sqrt(sum((Zt[i, .] - Zt[j, .]) :^ 2)))) {
		A[i, j] = 1
		A[j, i] = 1
	}
}
end
lsmt_mknet net40, edges("$LSM_REPO/dev/ref/random40_d2_edges.csv") n(40) dir(0)
lsmt_same chains, net(net40) opts(dim(2) seed(5) nsample(300) burnin(500) chains(3)) t1(1) t2(2)
lsmt_same chains5, net(net40) opts(dim(2) seed(5) nsample(300) burnin(500) chains(3)) t1(1) t2(5)
lsmt_same mklstarts, net(net40) opts(dim(2) seed(6) nsample(200) burnin(300) chains(2) mklstarts(4)) t1(1) t2(4)
lsmt_mknet net45, edges("$LSM_REPO/dev/ref/clu45_edges.csv") n(45) dir(0)
lsmt_same clusters, net(net45) opts(dim(2) seed(7) nsample(300) burnin(500) chains(3) clusters(3)) t1(1) t2(3)
* a larger network: one chain, the per-dyad passes and the MKL objective run on several threads
clear
set obs 200
set seed 3
forvalues j = 1/200 {
	quietly gen double __a`j' = .
}
mata: n = 200; Z = rnormal(n, 2, 0, 1); eta = -0.3 :+ rnormal(n, 1, 0, .6) :+ rnormal(1, n, 0, .6) :- lsmt_pd(Z); A = (runiform(n, n) :< invlogit(eta)); _diag(A, 0); st_store(., ., A); st_matrix("XC", rnormal(n, n, 0, 1))
nwset __a*, name(big) directed
lsmt_same inner, net(big) opts(dim(2) seed(8) nsample(40) burnin(100) pilot(2) plen(50) sender receiver edgecov(XC)) t1(1) t2(4)
* the persistent worker team of the per-dyad passes (networks of about 40 actors and more): n = 100, threads 1, 2, 4 and 8, and the same numbers as the version
* without a team (0.6.0: the fixed seed 7 on this network gave the checksum below)
clear
set obs 100
set seed 12345
forvalues j = 1/100 {
	quietly gen double a`j' = .
}
mata: n = 100; Zt = rnormal(n, 2, 0, 1.5); A = J(n, n, 0); for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) lsmt_tie(A, Zt, i, j); st_store(., ., A)
nwset a*, name(pool100) undirected
scalar ref = .
foreach th in 1 2 4 8 {
	quietly nwlsm pool100, dim(2) seed(7) nsample(1000) burnin(2000) threads(`th')
	mata: st_numscalar("h", sum(st_matrix("e(Z)") :* (1::rows(st_matrix("e(Z)")))) + 1000 * sum(st_matrix("e(b)")) + 1e6 * sum(st_matrix("e(summary)")))
	di as text "pool100 threads `th': checksum " %20.9f scalar(h)
	if scalar(ref) == . scalar ref = scalar(h)
	assert scalar(h) == scalar(ref)
}
assert abs(scalar(ref) - 3710932.303690526) < 1e-6
di as result "PASS pool100"
* errors
capture nwlsm big, threads(-1)
assert _rc == 198
capture nwlsm big, chains(0)
assert _rc == 198
capture nwlsm big, chains(2) nsample(10)
assert _rc == 198
capture nwlsm big, clusters(1)
assert _rc == 198
capture nwlsm big, clusters(7)
assert _rc == 198
display "PASS test_nwlsm_threads: identical results for every number of threads"
