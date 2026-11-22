* The default threads(0) (4 threads for one chain, 2 per chain for several chains, at most one per processor core) and threads(1), 2, 4 must give
* bit-identical e(b) and e(mkl_obj) for a given seed, with one chain and with four chains, on two directed 71-actor networks (Lazega).
do cscripts/nwlsm_testlib.do
capture nwwebuse lazega, nwclear
if _rc {
	display as text "lazega not loaded (rc `=_rc'), test skipped"
	exit
}
local worst = 0
foreach net in lazega_fr lazega_adv {
	foreach C in 1 4 {
		quietly nwlsm `net', dim(2) chains(`C') nsample(400) burnin(300) pilot(2) plen(50) seed(11) threads(1)
		matrix bref = e(b)
		local mref = e(mkl_obj)
		foreach th in 0 2 4 {
			if `th' == 0 quietly nwlsm `net', dim(2) chains(`C') nsample(400) burnin(300) pilot(2) plen(50) seed(11)
			else quietly nwlsm `net', dim(2) chains(`C') nsample(400) burnin(300) pilot(2) plen(50) seed(11) threads(`th')
			mata: st_local("db", strofreal(max(abs(st_matrix("bref") - st_matrix("e(b)")))))
			local dm = abs(`mref' - e(mkl_obj))
			di as text "`net' chains(`C') threads(`th') against threads(1): max|diff| b " `db' "  mkl_obj " `dm'
			assert `db' == 0 & `dm' == 0
		}
	}
}
di as result "PASS threads default"
