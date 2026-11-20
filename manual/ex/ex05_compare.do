* Comparing latent space models of the Lazega friendship network
clear all
set linesize 90
set scheme s1network
nwwebuse lazega, nwclear
log using out/compare_table.txt, text replace
quietly nwlsm lazega_fr, dim(2) seed(12)
estimates store d2
quietly nwlsm lazega_fr, dim(3) seed(12)
estimates store d3
quietly nwlsm lazega_fr, dim(2) absdiff(seniority) seed(12)
estimates store d2sen
quietly nwlsm lazega_fr, dim(2) sender receiver seed(12)
estimates store d2sr
nwcompare d2 d3 d2sen d2sr
log close
log using out/compare_select.txt, text replace
nwlsmselect lazega_fr, clusters(0 2 3 4) absdiff(seniority) seed(12)
log close
