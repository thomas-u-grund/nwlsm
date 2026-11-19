* Influence model: results, post-estimation, threads and chains, combinations and error paths.
do cscripts/nwlsm_testlib.do
do cscripts/nwlsm_infl_sim.do
lsmi_sim, n(40) dir(1) d(2) a0(1.5) b0(0.5) b1(0.6) b2(1) b3(-0.4) s2(0.25) g(8) wt(0) seed(31) zsd(2.5)
nwset __a*, name(ip) directed replace
* same seed gives the same fit; results do not depend on threads(); two chains
quietly nwlsm ip, dim(2) influence(y ylag x) neighbors(8) chains(2) threads(1) seed(9) nsample(400) burnin(500) interval(3)
matrix b1 = e(b)
matrix W1 = e(W)
matrix V1 = e(V)
assert e(chains) == 2 & "`e(influence)'" == "yes"
assert colsof(b1) == 5 & rowsof(e(rhat)) == 5 + 1 + 1 + 1
quietly nwlsm ip, dim(2) influence(y ylag x) neighbors(8) chains(2) threads(2) seed(9) nsample(400) burnin(500) interval(3)
assert mreldif(b1, e(b)) == 0 & mreldif(W1, e(W)) == 0 & mreldif(V1, e(V)) == 0
di as result "PASS threads(1) and threads(2) give identical results, chains(2)"
* structure of the results
matrix W = e(W)
mata: st_numscalar("__mx", max(abs(rowsum(st_matrix("W")) :- 1)))
assert scalar(__mx) < 1e-12
mata: st_numscalar("__nz", max(rowsum(st_matrix("W") :> 1e-12)))
* the weight matrix has the nonzero entries of draws: more than 8 positive entries per row in the average, none on the diagonal
mata: st_numscalar("__dg", max(abs(diagonal(st_matrix("W")))))
assert scalar(__dg) == 0
matrix PR = e(pred)
assert rowsof(PR) == 40 & colsof(PR) == 4
matrix cn = e(b)
assert "`: coleq cn'" == "network outcome outcome outcome outcome"
assert "`e(outcome)'" == "y" & "`e(lagvar)'" == "ylag" & "`e(weight)'" == "invdist" & e(neighbors) == 8
di as result "PASS weights rows sum to 1, no self weight, stored results"
* estat influence
estat influence, top(3) predict(pr) exposure(ex) exerted(xr)
assert reldif(r(influence), e(b)[1, 5]) < 1e-12
assert r(prob_positive) >= 0 & r(prob_positive) <= 1 & r(r2) > 0 & r(r2) <= 1 & r(r2) >= r(r2_noinfluence) - 1e-12
quietly summarize xr
assert reldif(r(sum), 40) < 1e-10
quietly correlate pr y
assert r(rho) > 0.5
estat influence, nolist graph
estat mcmcdiag, nolist
assert rowsof(r(diag)) == 5 + 1 + 1 + 1
di as result "PASS estat influence, estat mcmcdiag"
* the posterior mean coefficients agree with least squares on the posterior mean weights (loose: the weights are random, the prior vague)
matrix T = e(summary)
quietly gen double wy = .
mata: st_store(., "wy", st_matrix("W") * st_data(., "ylag"))
quietly regress y ylag x wy
matrix B = e(b)
assert abs(B[1, 4] - T[5, 1]) < 0.4
di as result "PASS posterior mean of the influence effect close to least squares on the posterior mean weights: " %6.3f T[5, 1] " vs " %6.3f B[1, 4]
* savedraws
tempfile dr
quietly nwlsm ip, dim(2) influence(y ylag x) neighbors(8) seed(9) nsample(200) burnin(300) interval(3) savedraws(`dr')
matrix bb = e(b)
preserve
use `dr', clear
confirm variable out_influence out_cons out_ylag out_x sigma2_y loglik_outcome beta_cons
quietly summarize out_influence
assert reldif(r(mean), bb[1, 5]) < 1e-9
restore
di as result "PASS savedraws"
* combinations: undirected + sociality + clusters + weights expdist; valued network; absdiff and edgecov
lsmi_sim, n(30) dir(0) d(2) a0(1.5) b0(0.5) b1(0.6) b2(1) b3(-0.4) s2(0.25) g(0) wt(1) seed(32) zsd(2.5)
nwset __a*, name(iu) undirected replace
quietly nwlsm iu, dim(2) sociality clusters(2) influence(y ylag x) weight(expdist) seed(3) nsample(300) burnin(400) interval(3) absdiff(x)
assert e(clusters) == 2 & e(neighbors) == 29 & "`e(weight)'" == "expdist"
estat influence, nolist
nwlsmplot, size(influence) name(lsmi_p)
nwlsmplot, size(influence) clusters name(lsmi_p2)
quietly nwlsm iu, dim(1) influence(y ylag) seed(3) nsample(200) burnin(200) interval(2)
mata: __tmpY = st_data(., ("y"))
nwset __a*, name(ivv) undirected valued replace
quietly nwlsm ivv, family(poisson) dim(2) influence(y ylag x) neighbors(5) seed(3) nsample(200) burnin(200) interval(2) obvar(25) osigma(2 1)
assert e(osigma_a) == 2 & e(obvar) == 25
di as result "PASS combinations: sociality, clusters, absdiff, expdist, one dimension, valued family"
* errors
nwset __a*, name(ie) undirected replace
foreach cmd in ///
	"nwlsm ie, influence(y ylag) bilinear" ///
	"nwlsm ie, influence(y)" ///
	"nwlsm ie, influence(y y)" ///
	"nwlsm ie, neighbors(3)" ///
	"nwlsm ie, influence(y ylag) weight(gauss)" ///
	"nwlsm ie, weight(expdist)" ///
	"nwlsm ie, influence(y ylag) obvar(-1)" ///
	"nwlsm ie, influence(y nosuchvar)" ///
	"nwlsm ie, influence(y ylag) neighbors(-2)" ///
	"nwlsm ie, influence(y ylag) osigma(1 -1)" {
	capture `cmd'
	assert _rc != 0
}
quietly replace ylag = . in 3
capture nwlsm ie, influence(y ylag)
assert _rc == 416
quietly replace ylag = 0.5 in 3
preserve
quietly drop in 1
capture nwlsm ie, influence(y ylag)
assert _rc == 198
restore
quietly nwlsm ie, dim(2) seed(1) nsample(100) burnin(100) interval(2)
capture estat influence
assert _rc == 198
capture nwlsmplot, size(influence)
assert _rc == 198
di as result "PASS error paths"
display "PASS test_nwlsm_infl_post"
