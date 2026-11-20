* The Lazega law firm: covariates, random effects, goodness of fit, diagnostics
clear all
set linesize 90
set scheme s1network
log using out/lazega_data.txt, text replace
nwwebuse lazega, nwclear
nwsummarize lazega_fr
describe, short
summarize seniority
log close
log using out/lazega_fit.txt, text replace
nwlsm lazega_fr, dim(2) absdiff(seniority) seed(12)
log close
log using out/lazega_select.txt, text replace
nwlsmselect lazega_fr, clusters(0 2 3 4) absdiff(seniority) seed(12)
log close
log using out/lazega_sr.txt, text replace
nwlsm lazega_fr, dim(2) sender receiver chains(2) seed(12)
log close
log using out/lazega_diag.txt, text replace
estat mcmcdiag
log close
log using out/lazega_gof.txt, text replace
estat gof, gof(distance) nsim(100) seed(1)
log close
log using out/lazega_edgecov.txt, text replace
nwtomata lazega_adv, mat(ADVM)
mata: _diag(ADVM, 0); st_matrix("ADV", ADVM)
nwlsm lazega_fr, dim(2) edgecov(ADV) seed(12)
log close
log using out/lazega_time.txt, text replace
timer clear
timer on 1
quietly nwlsm lazega_fr, dim(2) sender receiver seed(12) threads(1)
timer off 1
timer on 2
quietly nwlsm lazega_fr, dim(2) sender receiver seed(12) chains(4) threads(1)
timer off 2
timer on 3
quietly nwlsm lazega_fr, dim(2) sender receiver seed(12) chains(4) threads(4)
timer off 3
timer list
log close
* figures
quietly nwlsm lazega_fr, dim(2) sender receiver chains(2) seed(12)
nwlsmplot, color(office) colors(gray) size(sender) title("") xsize(5.5) ysize(5) name(l_map) export(../fig/fig_lazega_map.pdf)
estat gof, gof(idegree odegree distance) maxdeg(20) maxdist(8) nsim(100) seed(1) plot name(l_gof)
graph display l_gof, xsize(6.5) ysize(3.4)
graph export ../fig/fig_lazega_gof.pdf, name(l_gof) replace
