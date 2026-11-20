* The influence model on simulated data with a known influence effect
clear all
set linesize 90
set scheme s1network
log using out/infl_sim.txt, text replace
set obs 60
forvalues j = 1/60 {
    quietly generate double a`j' = .
}
quietly generate double y = .
quietly generate double ylag = .
quietly generate double x = .
mata:
n = 60; rseed(20261120)
Z = rnormal(n, 2, 0, 1.5)
D = J(n, n, 0)
for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) D[i, j] = sqrt(sum((Z[i, .] - Z[j, .]):^2))
A = (runiform(n, n) :< invlogit(1.5 :- D))
_diag(A, 0)
ylag = rnormal(n, 1, 0, 1)
x = rnormal(n, 1, 0, 1)
W = J(n, n, 0)
for (k = 1; k <= n; k++) {
    o = order(D[k, .]', 1)
    o = select(o, o :!= k)
    o = o[1..8]
    w = 1 :/ D[k, o]'
    W[k, o] = (w / sum(w))'
}
y = 0.5 :+ 0.6 * ylag :- 0.4 * x :+ 1 * (W * ylag) :+ rnormal(n, 1, 0, 0.5)
st_store(., (1..n), A)
st_store(., n + 1, y)
st_store(., n + 2, ylag)
st_store(., n + 3, x)
end
nwset a1-a60, name(sim) directed
nwsummarize sim
summarize y ylag x
log close
log using out/infl_fit.txt, text replace
nwlsm sim, dim(2) influence(y ylag x) neighbors(8) chains(2) nsample(2000) burnin(5000) seed(20261120)
log close
log using out/infl_estat.txt, text replace
estat influence, top(5) predict(pr) exposure(ex) exerted(xr)
log close
log using out/infl_diag.txt, text replace
estat mcmcdiag
log close
quietly nwlsm sim, dim(2) influence(y ylag x) neighbors(8) chains(2) nsample(2000) burnin(5000) seed(20261120)
nwlsmplot, size(influence) noedges title("") xsize(5.5) ysize(4.6) name(i_map) export(../fig/fig_infl_map.pdf)
estat influence, nolist graph name(i_pred)
graph display i_pred, xsize(5) ysize(4)
graph export ../fig/fig_infl_pred.pdf, name(i_pred) replace
