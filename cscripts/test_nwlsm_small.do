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

* degenerate networks (no tie / every tie): warning printed (red text in the log: checked through e(degenerate)), rc 0, finite e(b); an ordinary network: no warning
foreach dir in 1 0 {
	foreach full in 0 1 {
		clear
		set obs 8
		forvalues j = 1/8 {
			gen double __d`j' = `full'
			qui replace __d`j' = 0 in `j'
		}
		local dtype = cond(`dir', "directed", "undirected")
		nwset __d*, name(dg`dir'`full') `dtype'
		capture noisily nwlsm dg`dir'`full', dim(2) clusters(0) nsample(30) burnin(30) pilot(1) plen(10) seed(1)
		assert _rc == 0
		assert e(degenerate) == 1
		matrix __b = e(b)
		forvalues k = 1/`=colsof(__b)' {
			assert __b[1, `k'] < .
		}
		nwlsm
		nwdrop dg`dir'`full'
	}
}
clear
set obs 8
set seed 5
forvalues j = 1/8 {
	gen double __o`j' = runiform() < 0.4
	qui replace __o`j' = 0 in `j'
}
qui replace __o2 = 1 in 1
qui replace __o1 = 0 in 2
nwset __o*, name(dgord) directed
qui nwlsm dgord, dim(2) clusters(0) nsample(30) burnin(30) pilot(1) plen(10) seed(1)
assert e(degenerate) == 0
display as result "PASS degenerate networks"
