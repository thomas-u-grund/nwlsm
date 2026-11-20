* Sampson's monks: choosing the number of clusters, clusters against factions, valued ties
clear all
set linesize 90
set scheme s1network
log using out/monks_data.txt, text replace
nwwebuse sampson, nwclear
nwsummarize samplk3
nwgen choice3 = samplk3 > 0
nwsummarize choice3
log close
log using out/monks_bic.txt, text replace
nwlsm choice3, dim(2) seed(12)
estat bic
log close
log using out/monks_select.txt, text replace
nwlsmselect choice3, clusters(0 2 3 4) seed(12) refit graph gname(m_bic)
matrix list r(table), format(%9.2f)
display "best number of clusters: " r(best)
log close
log using out/monks_clusters.txt, text replace
estat clusters, generate(grp) nolist
tabulate grp_class faction
summarize grp_p1 grp_p2 grp_p3, separator(0)
generate pmax = max(grp_p1, grp_p2, grp_p3)
count if pmax < 0.8
log close
log using out/monks_effects.txt, text replace
nwlsmselect choice3, clusters(0 2 3 4) sender receiver seed(12)
nwlsm choice3, dim(2) clusters(3) sender receiver seed(12)
estat bic
display "starting points of the minimum Kullback-Leibler fit: " e(mklstarts)
nwlsm choice3, dim(2) clusters(3) sender receiver seed(12) mklstarts(1)
estat bic
estat bic, effobs(dyads)
estat bic, effobs(actors)
log close
log using out/monks_gof.txt, text replace
nwlsm choice3, dim(2) clusters(3) seed(12)
estat gof, gof(odegree distance) nsim(100) seed(1)
log close
log using out/monks_valued.txt, text replace
nwsummarize samplike
nwlsm samplike, family(binomial) trials(3) dim(2) clusters(3) seed(12)
estat gof, gof(distance values) nsim(100) seed(1)
estat gof, gof(distance) nsim(100) seed(1) cutoff(1)
estat ppcheck
log close
* figures
quietly nwlsm choice3, dim(2) clusters(3) seed(12)
nwlsmplot, clusters label(_nwnode) colors(gray) title("") xsize(5.5) ysize(4.6) name(m_clusters) export(../fig/fig_monks_clusters.pdf)
nwlsmplot, clusters clustermeans circles certainty label(_nwnode) colors(gray) title("") xsize(5.5) ysize(4.6) name(m_circles) export(../fig/fig_monks_circles.pdf)
nwlsmplot, color(faction) colors(gray) label(_nwnode) title("") xsize(5.5) ysize(4.6) name(m_faction) export(../fig/fig_monks_faction.pdf)
quietly nwlsmselect choice3, clusters(0 2 3 4) seed(12) graph gname(m_bic) gopts(title("") xsize(5.5) ysize(3.8))
graph export ../fig/fig_monks_bic.pdf, name(m_bic) replace
quietly nwlsm choice3, dim(2) clusters(3) seed(12)
estat ppcheck, graph name(m_ppcheck)
graph display m_ppcheck, xsize(5) ysize(4)
graph export ../fig/fig_monks_ppcheck.pdf, name(m_ppcheck) replace
