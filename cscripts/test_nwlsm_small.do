* nwlsm on very small networks (3 to 6 actors): the optimiser's work space once overflowed when the number of position
* coordinates was below 10. Must run without a crash and give finite results, directed and undirected, several seeds.
do cscripts/nwlsm_testlib.do
local nruns 0
foreach n in 3 4 5 6 {
	foreach dir in 1 0 {
		foreach s in 1 2 3 {
			set seed `=100 * `n' + 10 * `dir' + `s''
			clear
			set obs `n'
			forvalues j = 1/`n' {
				gen double __a`j' = runiform() < 0.5
				qui replace __a`j' = 0 in `j'
			}
			* make sure the network has at least one tie
			qui replace __a`=mod(`s', `n') + 1' = 1 in `=mod(`s' + 1, `n') + 1'
			local dtype = cond(`dir', "directed", "undirected")
			nwset __a*, name(sm`n'`dir'`s') `dtype'
			qui nwlsm sm`n'`dir'`s', dim(2) clusters(0) nsample(30) burnin(30) pilot(1) plen(10) seed(`s')
			assert e(mkl_obj) < .
			assert e(mkl_iter) >= 0 & e(mkl_iter) < .
			matrix __b = e(b)
			forvalues k = 1/`=colsof(__b)' {
				assert __b[1, `k'] < .
			}
			local ++nruns
			nwdrop sm`n'`dir'`s'
		}
	}
}
assert `nruns' == 24
di "PASS test_nwlsm_small (`nruns' runs)"
