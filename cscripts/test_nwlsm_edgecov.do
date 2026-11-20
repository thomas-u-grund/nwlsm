* edgecov(): the diagonal of a one-mode network is not a dyad. A network (or a matrix from nwtomata) has a missing diagonal; the fit must equal
* the fit with the diagonal set to 0; missing off-diagonal values stay an error (r(416)); the matrix is read as before for directed and undirected networks.
do cscripts/nwlsm_testlib.do
set seed 21
clear
set obs 14
forvalues j = 1/14 {
	gen double __a`j' = runiform() < 0.3
	gen double __c`j' = runiform() < 0.35
}
local fast "dim(2) nsample(300) burnin(300) pilot(2) plen(50) seed(7)"
foreach dir in directed undirected {
	preserve
	nwset __a*, name(y_`dir') `dir'
	nwset __c*, name(x_`dir') `dir'
	restore
	* the diagonal of the network matrix is missing
	nwtomata x_`dir', mat(XM)
	mata: assert(XM[1, 1] >= .)
	mata: X0 = XM; for (i = 1; i <= 14; i++) X0[i, i] = 0; st_matrix("XM0", X0)
	mata: st_matrix("XMm", XM)
	quietly nwlsm y_`dir', edgecov(x_`dir') `fast'
	matrix b_net = e(b)
	quietly nwlsm y_`dir', edgecov(XMm) `fast'
	assert mreldif(b_net, e(b)) == 0
	quietly nwlsm y_`dir', edgecov(XM0) `fast'
	assert mreldif(b_net, e(b)) == 0
	* an off-diagonal missing value is still an error
	mata: X1 = X0; X1[2, 3] = .; st_matrix("XM1", X1)
	capture nwlsm y_`dir', edgecov(XM1) `fast'
	assert _rc == 416
}
* the Lazega lawyers (network adviser, diagonal missing), if the data can be loaded
capture nwwebuse lazega, nwclear
if _rc == 0 {
	capture nwlsm lazega_fr, dim(2) edgecov(lazega_adv) nsample(200) burnin(200) pilot(2) plen(50) seed(3)
	assert _rc == 0
}
else display "lazega not loaded (rc `=_rc'), network part skipped"
display "PASS test_nwlsm_edgecov"
