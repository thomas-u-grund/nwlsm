* Zachary's karate club (34 members, 78 ties), loaded from the nwcommands website
clear all
set linesize 90
set scheme s1network
log using out/version.txt, text replace
nwlsm version
log close
log using out/karate_data.txt, text replace
nwwebuse zachary, nwclear
nwsummarize zachary_bin
log close
log using out/karate_fit.txt, text replace
nwlsm zachary_bin, dim(2) nsample(2000) burnin(5000) interval(10) chains(2) seed(20261120)
log close
log using out/karate_diag.txt, text replace
estat mcmcdiag
log close
log using out/karate_pos.txt, text replace
estat positions, generate(pos) nolist
list _nwnode pos1 pos2 in 1/8, noobs
log close
log using out/karate_gof.txt, text replace
estat gof, gof(degree distance) nsim(100) seed(1)
log close
log using out/karate_ppcheck.txt, text replace
estat ppcheck
log close
log using out/karate_replay.txt, text replace
nwlsm
display "version: " e(version)
log close
* the same model with a random sociality effect and with the bilinear term
log using out/karate_soc.txt, text replace
nwlsm zachary_bin, dim(2) sociality nsample(2000) burnin(5000) chains(2) seed(20261120)
log close
log using out/karate_bil.txt, text replace
nwlsm zachary_bin, dim(2) bilinear nsample(2000) burnin(5000) chains(2) seed(20261120)
log close
* threads: the same seed gives the same numbers for any number of threads
log using out/karate_threads.txt, text replace
quietly nwlsm zachary_bin, dim(2) chains(2) threads(1) seed(5) nsample(1000) burnin(2000)
matrix b1 = e(b)
display "threads(1): _cons = " %9.6f b1[1,1]
quietly nwlsm zachary_bin, dim(2) chains(2) threads(2) seed(5) nsample(1000) burnin(2000)
matrix b2 = e(b)
display "threads(2): _cons = " %9.6f b2[1,1]
log close
* saved draws and the draws in memory
log using out/karate_draws.txt, text replace
quietly nwlsm zachary_bin, dim(2) seed(5) nsample(500) burnin(2000) savedraws(draws_karate.dta)
preserve
use draws_karate.dta, clear
describe, short
summarize beta_cons loglik, separator(0)
restore
mata: st_numscalar("nd", rows(LsmSample("beta")))
display "draws of beta in Mata: " nd
log close
erase draws_karate.dta
* figures: map, map with draws, trace plots, goodness of fit, calibration
quietly nwlsm zachary_bin, dim(2) nsample(2000) burnin(5000) interval(10) chains(2) seed(20261120)
nwlsmplot, title("") xsize(5.5) ysize(4.6) name(k_map) export(../fig/fig_karate_map.pdf)
nwlsmplot, draws(100) title("") xsize(5.5) ysize(4.6) name(k_draws) export(../fig/fig_karate_draws.pdf)
estat mcmcdiag, graph nolist name(k_trace)
graph display k_trace, xsize(6.5) ysize(6.5)
graph export ../fig/fig_karate_trace.pdf, name(k_trace) replace
