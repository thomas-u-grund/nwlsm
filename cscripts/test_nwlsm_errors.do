* Error paths and e() results of nwlsm.
do cscripts/nwlsm_testlib.do
set seed 11
clear
set obs 12
forvalues j = 1/12 {
	gen double __a`j' = runiform() < 0.3
}
nwset __a*, name(g) directed
mata: st_matrix("__BADM", J(5, 5, 1)); st_matrix("__OKM", abs(runiform(12, 1) * J(1, 12, 1) - J(12, 1, 1) * runiform(12, 1)'))
local fast "nsample(50) burnin(100) pilot(1) plen(20)"
* dim
capture nwlsm g, dim(0) `fast'
assert _rc == 198
capture nwlsm g, dim(11) `fast'
assert _rc == 198
capture nwlsm g, dim(2.5) `fast'
assert _rc != 0
* sampler settings
capture nwlsm g, nsample(5) `fast'
assert _rc == 198
capture nwlsm g, interval(0) nsample(50) burnin(10)
assert _rc == 198
capture nwlsm g, burnin(-1) nsample(50)
assert _rc == 198
capture nwlsm g, seed(-3) `fast'
assert _rc == 198
capture nwlsm g, seed(abc) `fast'
assert _rc == 198
* not a network
capture nwlsm nosuchnet, `fast'
assert _rc != 0
capture nwlsm, `fast'
assert _rc != 0
* covariates
capture nwlsm g, edgecov(__BADM) `fast'
assert _rc == 198
capture nwlsm g, edgecov(nosuchthing) `fast'
assert _rc == 198
capture nwlsm g, absdiff(nosuchvar) `fast'
assert _rc != 0
capture nwlsm g, absdiff(__a1) `fast'
assert _rc == 0
capture nwlsm g, edgecov(__OKM __OKM) `fast'
assert _rc == 198
* replay without estimates
ereturn clear
capture nwlsm
assert _rc == 301
* a good fit: e() results and replay
quietly nwlsm g, edgecov(__OKM) dim(2) seed(5) nsample(200) burnin(200) pilot(2) plen(50)
assert "`e(cmd)'" == "nwlsm"
assert e(dim) == 2 & e(nsample) == 200 & e(N_nodes) == 12 & e(seed) == 5
assert colsof(e(b)) == 2 & rowsof(e(V)) == 2 & rowsof(e(D)) == 12 & rowsof(e(Z)) == 12 & colsof(e(Z)) == 2
assert "`: colnames e(b)'" == "_cons edgecov___OKM"
assert abs(e(D)[1, 2] - e(D)[2, 1]) < 1e-12 & e(D)[3, 3] == 0
assert e(sigma2_mean) > 0 & e(prior_df) > 0
nwlsm
* same seed, same result
matrix b1 = e(b)
quietly nwlsm g, edgecov(__OKM) dim(2) seed(5) nsample(200) burnin(200) pilot(2) plen(50)
assert mreldif(b1, e(b)) == 0
display "PASS test_nwlsm_errors"
