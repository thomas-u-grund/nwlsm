* estat mcmcdiag: (1) the Geweke flag uses the Bonferroni-corrected normal quantile for parameters x chains, so a converged run with 4 or 8 chains
* is not flagged (the largest |z| over many chains and parameters is routinely above 2); a short run is flagged (ESS or Geweke);
* (2) no output of estat mcmcdiag, nwlsm, estat gof/ppcheck/bic/clusters or nwlsmselect wraps at linesize 80.
do cscripts/nwlsm_testlib.do
global REF "$LSM_REPO/dev/ref"
lsmt_mknet und, edges("$REF/clu45_edges.csv") n(45) dir(0)
lsmt_mknet dir, edges("$REF/sampson_d2_edges.csv") n(18) dir(1)

* (1) converged runs, 4 and 8 chains, with several parameters
foreach C in 1 4 8 {
	quietly nwlsm und, dim(2) sociality chains(`C') nsample(4000) burnin(20000) seed(5)
	estat mcmcdiag
	local np = rowsof(r(diag))
	local k = `np' * `C'
	assert reldif(r(geweke_crit), invnormal(1 - 0.05 / (2 * `k'))) < 1e-12
	assert r(geweke_ntests) == `k'
	assert r(geweke_max) <= r(geweke_crit)
	assert r(geweke_flag) == 0
	assert r(ess_flag) == 0
	* the maximum is the maximum of the printed column
	mata: st_numscalar("__m", max(abs(st_matrix("r(diag)")[., 2])))
	assert reldif(r(geweke_max), scalar(__m)) < 1e-12
	di "PASS: converged run, `C' chains: max |z| " %5.2f r(geweke_max) " <= " %5.2f r(geweke_crit)
}
* the old rule |z| > 2 would have flagged a row of this 8-chain run or is at least not the rule any more
assert r(geweke_crit) > 2

* a deliberately short run is flagged
quietly nwlsm dir, dim(2) nsample(60) burnin(0) interval(1) seed(3)
estat mcmcdiag
assert r(geweke_flag) == 1 | r(ess_flag) == 1
assert r(ess_min) < 100 | r(geweke_flag) == 1
di "PASS: short run flagged (geweke_flag " r(geweke_flag) ", ess_flag " r(ess_flag) ")"

* (2) no wrapped output at linesize 80
capture log close __wrap
log using "$LSM_SCRATCH/wrap_mcmcdiag.log", text replace name(__wrap)
set linesize 80
quietly nwlsm dir, dim(2) sender receiver sociality chains(2) nsample(300) burnin(300) seed(2)
estat mcmcdiag
estat gof, nsim(10) seed(1)
estat ppcheck
estat bic
quietly nwlsm und, dim(2) clusters(2) nsample(300) burnin(300) seed(2)
estat clusters
estat mcmcdiag
nwlsmselect und, clusters(0 2) nsample(200) burnin(200) seed(2)
log close __wrap
set linesize 255
tempname fh
file open `fh' using "$LSM_SCRATCH/wrap_mcmcdiag.log", read text
local incmd 0
local bad 0
local width 0
file read `fh' line
while r(eof) == 0 {
	if substr(`"`line'"', 1, 2) == ". " local incmd 1
	else if substr(`"`line'"', 1, 2) == "> " {
		if !`incmd' & strpos(`"`line'"', "mcdiag.log") == 0 {  // the log header and footer print the path of the log file
			local ++bad
			di as error "wrapped line: `line'"
		}
	}
	else local incmd 0
	file read `fh' line
}
file close `fh'
assert `bad' == 0
di "PASS: no wrapped output at linesize 80"
