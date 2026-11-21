* estat mcmcdiag: (1) lsm_geweke() reproduces coda::geweke.diag() of R (values computed with coda 0.19-4.1 for fixed series); (2) the printed
* Geweke z, r(geweke_z) (chains x parameters), the count of |z| > 1.96, the nominal flag; the number exceeding 1.96 in converged runs is within
* a binomial plausibility band; a clearly non-converged run is flagged; (3) no output of estat mcmcdiag, nwlsm, estat gof/ppcheck/bic/clusters
* or nwlsmselect wraps at linesize 80.
do cscripts/nwlsm_testlib.do
global REF "$LSM_REPO/dev/ref"
lsmt_mknet und, edges("$REF/clu45_edges.csv") n(45) dir(0)
lsmt_mknet dir, edges("$REF/sampson_d2_edges.csv") n(18) dir(1)

* (1) numerical check against coda::geweke.diag(x, frac1 = 0.1, frac2 = 0.5). Series: x_t = rho x_{t-1} + h_t + slope t, with the
* deterministic noise h_t = frac(sin(12.9898 t) 43758.5453) - 0.5. R: geweke.diag(mcmc(x))$z
mata:
real colvector lsmt_series(real scalar S, real scalar rho, real scalar slope)
{
	real colvector t, h, x
	real scalar k
	t = (1::S)
	h = sin(t * 12.9898) * 43758.5453
	h = h - floor(h) :- 0.5
	x = h
	for (k = 2; k <= S; k++) x[k] = rho * x[k - 1] + h[k]
	return(x + slope * t)
}
end
mata: st_numscalar("__z1", lsm_geweke(lsmt_series(800, 0.8, 0.002)))
mata: st_numscalar("__z2", lsm_geweke(lsmt_series(1000, 0, 0)))
mata: st_numscalar("__z3", lsm_geweke(lsmt_series(333, 0.95, 0)))
assert abs(scalar(__z1) - (-5.53163191535)) < 1e-6
assert abs(scalar(__z2) - (-1.04119064273)) < 1e-6
assert abs(scalar(__z3) - (-0.495064410572)) < 1e-6
di "PASS: lsm_geweke() equals coda::geweke.diag() (" %9.6f scalar(__z1) ", " %9.6f scalar(__z2) ", " %9.6f scalar(__z3) ")"
* a series with a perfect linear trend has zero spectral density in coda (z missing); a constant too
mata: st_numscalar("__z4", lsm_geweke((1::200) * 0.5))
mata: st_numscalar("__z5", lsm_geweke(J(200, 1, 3)))
assert scalar(__z4) >= . & scalar(__z5) >= .
di "PASS: degenerate series give a missing z"

* (2) converged runs, 1, 4 and 8 chains: counts and semantics
local totk 0
local totex 0
foreach C in 1 4 8 {
	quietly nwlsm und, dim(2) sociality chains(`C') nsample(4000) burnin(20000) seed(5)
	estat mcmcdiag
	local np = rowsof(r(diag))
	local k = `np' * `C'
	assert r(geweke_ntests) == `k'
	* r(geweke_z): chains x parameters; the count and the maximum come from it; the printed column is the value of largest size per row
	matrix Z = r(geweke_z)
	assert rowsof(Z) == `C' & colsof(Z) == `np'
	mata: Z = st_matrix("Z"); st_numscalar("__ex", sum(abs(Z) :> 1.959963984540054)); st_numscalar("__mx", max(abs(Z)))
	assert r(geweke_nexceed) == scalar(__ex)
	assert reldif(r(geweke_max), scalar(__mx)) < 1e-12
	mata: D = st_matrix("r(diag)"); st_numscalar("__ok", all(abs(D[., 2]) :== rowmax(abs(Z'))))
	assert scalar(__ok) == 1
	* plausibility. One chain of 4000 draws: the count is within the binomial band of 5% (upper tail probability above 0.1%). Several chains
	* are shorter (500 draws with 8 chains: the first 10% are 50 draws, whose autoregressive spectral density estimate is noisy), so in
	* settled runs more than 5% of the z values exceed 1.96 (16 to 26% for 4 and 8 chains in 48 settled fits): less than half and none extreme.
	if `C' == 1 assert binomialtail(`k', r(geweke_nexceed), 0.05) > 0.001
	else assert r(geweke_nexceed) < `k' / 2
	assert r(geweke_max) < 6
	assert r(ess_flag) == 0
	local totk = `totk' + `k'
	local totex = `totex' + r(geweke_nexceed)
	di "PASS: converged run, `C' chains: " r(geweke_nexceed) " of `k' |z| exceed 1.96 (about " %3.1f 0.05 * `k' " expected), max |z| " %5.2f r(geweke_max)
}
assert `totex' < `totk' / 2
di "PASS: " `totex' " of " `totk' " z values exceed 1.96 over the three runs (about " %3.1f 0.05 * `totk' " expected)"

* a deliberately short run: the ESS warning (the Geweke z is informational, there is no Geweke flag)
quietly nwlsm dir, dim(2) nsample(60) burnin(0) interval(1) seed(3)
estat mcmcdiag
assert r(ess_flag) == 1 & r(ess_min) < 100
capture confirm scalar r(geweke_flag)
assert _rc != 0
di "PASS: short run: ESS warning (ess_flag " r(ess_flag) "), no Geweke flag"
* a run that is far from converged (a single chain of 400 draws that has hardly left its start): all three z values far beyond 1.96
quietly nwlsm und, dim(2) chains(1) nsample(400) burnin(100) interval(1) seed(3)
estat mcmcdiag
assert r(geweke_nexceed) == 3 & r(geweke_max) > 5
di "PASS: non-converged run: " r(geweke_nexceed) " of " r(geweke_ntests) " |z| exceed 1.96, max " %5.2f r(geweke_max)

* (3) no wrapped output at linesize 80
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
