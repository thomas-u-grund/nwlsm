* Post-estimation and options of version 0.3.0: estat positions, estat mcmcdiag, estat predict, nwlsmplot, savedraws with the variance draws,
* the checks of the effect options, and the effective sample size / Geweke statistics on series with known properties.
do cscripts/nwlsm_testlib.do
set seed 31
clear
set obs 14
forvalues j = 1/14 {
	gen double __a`j' = runiform() < 0.3
}
mata: A = st_data(., .); A = uppertriangle(A, 0); A = A + A'; _diag(A, 0); st_store(., ., A)
nwset __a*, name(u) undirected
clear
set obs 14
forvalues j = 1/14 {
	gen double __a`j' = runiform() < 0.3
}
nwset __a*, name(g) directed
local fast "nsample(300) burnin(300) pilot(2) plen(50)"
* ---- error paths of the new options
capture nwlsm u, sender `fast'
assert _rc == 198
capture nwlsm u, receiver `fast'
assert _rc == 198
capture nwlsm u, sociality sender `fast'
assert _rc == 198
capture nwlsm g, sender sendervar(0) `fast'
assert _rc == 198
capture nwlsm g, sociality socialitydf(-1) `fast'
assert _rc == 198
capture nwlsm g, mklstarts(0) `fast'
assert _rc == 198
capture nwlsm g, mklstarts(2.5) `fast'
assert _rc != 0
capture estat positions
assert _rc == 301 | _rc == 321
* ---- a directed fit with sender and receiver
quietly nwlsm g, dim(2) sender receiver seed(3) `fast'
assert "`e(effects)'" == "sender receiver" & "`e(estat_cmd)'" == "nwlsm_estat"
assert rowsof(e(vsummary)) == 3 & rowsof(e(sender)) == 14 & colsof(e(sender)) == 2 & rowsof(e(receiver)) == 14
assert rowsof(e(P)) == 14 & rowsof(e(mkl_effects)) == 14 & colsof(e(mkl_b)) == 1
assert e(mkl_obj) < 0 & e(mkl_iter) > 0 & e(mkl_start) >= 1 & e(mkl_start) <= 300
* the minimum-KL effects are centred (the means are in the intercept)
mata: st_local("m1", strofreal(abs(mean(st_matrix("e(mkl_effects)")[., 1])) + abs(mean(st_matrix("e(mkl_effects)")[., 2]))))
assert `m1' < 1e-10
* probabilities are probabilities and symmetric information is not imposed in the directed case
mata: st_local("pmin", strofreal(min(st_matrix("e(P)")))); mata: st_local("pmax", strofreal(max(st_matrix("e(P)"))))
assert `pmin' >= 0 & `pmax' <= 1
estat positions
matrix Z1 = r(Z)
assert mreldif(Z1, e(Z)) == 0
estat positions, mean
matrix Z2 = r(Z)
assert mreldif(Z2, e(Zmean)) == 0
* the mean of the aligned draws and the minimum-KL positions describe the same configuration up to shrinkage: their distances are correlated
mata: st_local("cz", strofreal(correlation((vec(lsmt_pd(st_matrix("e(Z)"))), vec(lsmt_pd(st_matrix("e(Zmean)")))))[1, 2]))
assert `cz' > 0.7
estat mcmcdiag
matrix DG = r(diag)
assert rowsof(DG) == 5 & colsof(DG) == 5
forvalues k = 1/5 {
	assert DG[`k', 1] >= 1 & DG[`k', 1] <= 300 & abs(DG[`k', 2]) < 10
}
assert "`: rownames DG'" == "_cons loglik sigma2 sender_var receiver_var" | "`: rownames DG'" == "_cons loglik sigma2 sender_var receiver_var "
capture drop pos1 pos2
clear
set obs 14
estat positions, generate(pos) nolist
confirm variable pos1 pos2
assert abs(pos1[3] - el(e(Z), 3, 1)) < 1e-12
capture estat positions, generate(pos) nolist
assert _rc == 110
set obs 15
capture estat positions, generate(zz)
assert _rc == 198
clear
* predict: a valued network of posterior mean tie probabilities
estat predict, generate(pnet)
matrix PP = r(P)
assert mreldif(PP, e(P)) == 0
_nwsyntax pnet, max(1)
assert "`valued'" == "true" & "`directed'" == "true" & `nodes' == 14
mata: st_local("dd", strofreal(max(abs(*(`netobj'->get_matrix()) - st_matrix("e(P)"))) ))
assert `dd' < 1e-12
capture estat predict, generate(pnet)
assert _rc == 110
estat predict, generate(pnet) replace
capture estat nosuch
assert _rc == 321
* ---- savedraws holds the variance draws
tempfile dr
quietly nwlsm g, dim(2) sender receiver seed(3) `fast' savedraws(`dr')
preserve
use `dr', clear
confirm variable beta_cons loglik sigma2 sender_var receiver_var
assert _N == 300 & sender_var > 0 & receiver_var > 0
restore
* ---- an undirected fit with sociality; bilinear; one-dimensional
quietly nwlsm u, dim(2) sociality seed(4) `fast'
assert "`e(effects)'" == "sociality" & rowsof(e(sociality)) == 14 & rowsof(e(vsummary)) == 2
estat mcmcdiag, nolist
assert rowsof(r(diag)) == 4
quietly nwlsm u, dim(2) bilinear seed(4) `fast'
assert "`e(latent)'" == "bilinear"
estat positions, nolist
quietly nwlsm u, dim(1) seed(4) `fast' sociality
estat positions, nolist
assert colsof(r(Z)) == 1
capture nwlsmplot
assert _rc == 198
* ---- the plot
quietly nwlsm g, dim(2) sender receiver seed(3) `fast'
local out "$LSM_SCRATCH"
nwlsmplot, name(lsm_p1)
graph export "`out'/lsm_p1.png", replace width(600)
confirm file "`out'/lsm_p1.png"
clear
set obs 14
gen byte grp = mod(_n, 3)
gen str6 nm = "a" + string(_n)
gen double sz = _n
nwlsmplot, color(grp) size(sender) label(nm) name(lsm_p2) dims(1 2) title("Test")
graph export "`out'/lsm_p2.png", replace width(600)
confirm file "`out'/lsm_p2.png"
nwlsmplot, color(grp) size(sz) noedges mean name(lsm_p3)
capture nwlsmplot, size(sociality)
assert _rc == 198
capture nwlsmplot, size(nosuchvar)
assert _rc == 198
capture nwlsmplot, dims(1 3)
assert _rc == 198
capture nwlsmplot, dims(1 1)
assert _rc == 198
set obs 20
capture nwlsmplot, color(grp)
assert _rc == 198
* ---- effective sample size and Geweke statistic on series with known properties
mata:
mata set matastrict off
void lsmt_series()
{
	real scalar S, rho, k, ess, z
	real colvector x, e

	rseed(77)
	S = 20000
	e = rnormal(S, 1, 0, 1)
	// independent draws: ESS close to S
	ess = lsm_ess(e)
	printf("iid: ESS %f of %f\n", ess, S)
	if (ess < 0.8 * S) exit(9)
	// AR(1) with rho = 0.8: ESS = S (1 - rho) / (1 + rho) = 0.111 S
	rho = 0.8
	x = e
	for (k = 2; k <= S; k++) x[k] = rho * x[k - 1] + e[k]
	ess = lsm_ess(x)
	printf("AR(1), rho 0.8: ESS %f, theory %f\n", ess, S * (1 - rho) / (1 + rho))
	if (abs(ess / (S * (1 - rho) / (1 + rho)) - 1) > 0.3) exit(9)
	// Geweke: no trend -> small, trend -> large
	z = lsm_geweke(e)
	printf("Geweke z, iid: %f\n", z)
	if (abs(z) > 4) exit(9)
	z = lsm_geweke(e + (1::S) / S * 3)
	printf("Geweke z, trend: %f\n", z)
	if (abs(z) < 5) exit(9)
}
lsmt_series()
end
display "PASS test_nwlsm_post: estat positions, mcmcdiag, predict, savedraws, the plot, option checks, ESS and Geweke statistics"
* chains: Gelman-Rubin statistics and the diagnostics table with several chains; estat clusters and the plot by cluster
quietly nwlsm g, dim(2) chains(3) nsample(300) burnin(500) seed(4)
assert e(chains) == 3 & e(nsample) == 300 & e(nsample_chain) == 100
matrix RH = e(rhat)
assert rowsof(RH) == 3 & colsof(RH) == 1
mata: assert(all(st_matrix("RH") :> 0.5) & all(st_matrix("RH") :< 3))
estat mcmcdiag
matrix DG = r(diag)
assert rowsof(DG) == 3 & colsof(DG) == 5 & DG[1, 5] < .
quietly nwlsm g, dim(2) clusters(2) nsample(300) burnin(500) seed(4)
estat clusters, nolist
matrix QQ = r(Q)
assert rowsof(QQ) == 14 & colsof(QQ) == 2
clear
set obs 14
estat clusters, generate(cl) nolist
assert cl_class < . & inlist(cl_class, 1, 2)
nwlsmplot, clusters name(gcl)
capture nwlsmplot, clusters color(cl_class)
assert _rc == 198
quietly nwlsm g, dim(2) nsample(300) burnin(500) seed(4)
capture estat clusters
assert _rc == 198
capture nwlsmplot, clusters
assert _rc == 198
display "PASS test_nwlsm_post: chains and clusters post-estimation"
