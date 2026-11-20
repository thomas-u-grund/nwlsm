* Plot options of nwlsmplot (arrows, clustermeans, circles, draws(), certainty, labsize/labposition, noaxes, export()) and the graph option of nwlsmselect,
* on an undirected fit (planted-cluster network, 45 actors) and a directed fit (Sampson monks, 18 actors); error paths; the data in memory are unchanged.
do cscripts/nwlsm_testlib.do
set linesize 255
global REF "$LSM_REPO/dev/ref"

capture program drop lsmt_gexists
program lsmt_gexists
	args nm
	capture graph describe `nm'
	if _rc {
		di as error "FAIL: graph `nm' does not exist"
		exit 9
	}
	graph drop `nm'
end

lsmt_mknet und, edges("$REF/clu45_edges.csv") n(45) dir(0)
lsmt_mknet dir, edges("$REF/sampson_d2_edges.csv") n(18) dir(1)

* ---- undirected fit with clusters
quietly nwlsm und, dim(2) clusters(2) seed(3) nsample(400) burnin(500)
quietly set obs 45
quietly gen str6 nm = "a" + string(_n)
quietly gen double grp = mod(_n, 3)
local N0 = _N
nwlsmplot, name(g1)
lsmt_gexists g1
nwlsmplot, clusters clustermeans circles name(g2)
lsmt_gexists g2
nwlsmplot, clusters clustermeans circles certainty draws(60) name(g3)
lsmt_gexists g3
nwlsmplot, clusters draws(30) noedges name(g4)
lsmt_gexists g4
nwlsmplot, draws(25) mean name(g5)
lsmt_gexists g5
nwlsmplot, color(grp) draws(20) label(nm) labsize(small) labposition(6) name(g6)
lsmt_gexists g6
nwlsmplot, noaxes clusters name(g7) title("no axes") scheme(s1mono)
lsmt_gexists g7
nwlsmplot, noarrows name(g8)
lsmt_gexists g8
nwlsmplot, clusters circles clustermeans export("$LSM_SCRATCH/test_nwlsm_plot.png", width(600)) name(g9)
lsmt_gexists g9
confirm file "$LSM_SCRATCH/test_nwlsm_plot.png"
nwlsmplot, saving("$LSM_SCRATCH/test_nwlsm_plot.gph", replace) name(g10)
lsmt_gexists g10
confirm file "$LSM_SCRATCH/test_nwlsm_plot.gph"
capture nwlsmplot, draws(5000)
assert _rc == 198
capture nwlsmplot, draws(-1)
assert _rc == 198
capture nwlsmplot, arrows
assert _rc == 198
capture nwlsmplot, certainty
assert _rc == 198
capture nwlsmplot, labsize(small)
assert _rc == 198
capture nwlsmplot, arrows noarrows
assert _rc == 198
assert _N == `N0'
describe, short

* ---- equal scale: aspect = y range / x range, equal padding, so one unit has the same length on both axes; noequal and aspectratio() work
quietly nwlsm und, dim(2) clusters(2) seed(3) nsample(400) burnin(500)
nwlsmplot, clusters circles clustermeans draws(30) name(e1)
assert abs(r(aspect) - (r(ymax) - r(ymin)) / (r(xmax) - r(xmin))) < 1e-9
assert r(xmax) > r(xmin) & r(ymax) > r(ymin)
* the range includes the circles: mean +- radius lies inside
matrix CM = e(cmeans)
matrix VS = e(vsummary)
forvalues g = 1/2 {
	assert CM[`g', 1] - sqrt(VS[`g', 1]) >= r(xmin) & CM[`g', 1] + sqrt(VS[`g', 1]) <= r(xmax)
	assert CM[`g', 2] - sqrt(VS[`g', 1]) >= r(ymin) & CM[`g', 2] + sqrt(VS[`g', 1]) <= r(ymax)
}
lsmt_gexists e1
nwlsmplot, noequal clusters circles name(e2)
assert r(aspect) == .
lsmt_gexists e2
nwlsmplot, aspectratio(0.5) name(e3)
lsmt_gexists e3
nwlsmplot, noaxes name(e4)
lsmt_gexists e4

* ---- fit without clusters: cluster options are errors
quietly nwlsm und, dim(2) seed(3) nsample(300) burnin(300)
foreach o in clustermeans circles "clusters certainty" {
	capture nwlsmplot, `o'
	assert _rc == 198
}
* ---- a later fit replaces the draws: draws() belong to the current fit
quietly nwlsm und, dim(3) seed(4) nsample(300) burnin(300)
nwlsmplot, dims(1 3) draws(20) name(g12)
lsmt_gexists g12
quietly nwlsm und, dim(2) seed(3) nsample(300) burnin(300)

* ---- draws not in memory: mismatch is detected (external from another fit size)
mata: mata drop __nwlsm_Z
capture nwlsmplot, draws(10)
assert _rc == 498
quietly nwlsm und, dim(2) seed(3) nsample(300) burnin(300)

* ---- directed fit: arrows by default
clear
quietly nwlsm dir, dim(2) seed(5) nsample(400) burnin(500)
nwlsmplot, name(d1)
lsmt_gexists d1
nwlsmplot, arrows draws(40) name(d2)
lsmt_gexists d2
nwlsmplot, noarrows name(d3)
lsmt_gexists d3
nwlsmplot, noedges draws(15) name(d4)
lsmt_gexists d4
capture nwlsmplot, arrows noedges
assert _rc == 198
quietly nwlsm dir, dim(2) clusters(2) seed(5) nsample(400) burnin(500)
nwlsmplot, clusters clustermeans circles certainty draws(50) name(d5)
lsmt_gexists d5
quietly set obs 18
quietly gen str6 nm = "m" + string(_n)
nwlsmplot, clusters circles arrows label(nm) labsize(tiny) name(d6)
lsmt_gexists d6
nwlsmplot, clusters label(nm) labposition(auto) noequal name(d6b)
lsmt_gexists d6b
nwlsmplot, clusters label(nm) labposition(3) name(d6c)
lsmt_gexists d6c
* undirected network: arrows are refused
quietly nwlsm und, dim(2) seed(3) nsample(300) burnin(300)
capture nwlsmplot, arrows
assert _rc == 198

* ---- nwlsmselect, graph
nwlsmselect und, clusters(0 2 3) seed(2) nsample(300) burnin(300) nolist graph gname(s1) gopts(title("BIC test"))
matrix T = r(table)
local bst = r(best)
lsmt_gexists s1
assert rowsof(T) == 3
nwlsmselect dir, clusters(0 2) seed(2) nsample(300) burnin(300) nolist graph sender receiver gname(s2)
lsmt_gexists s2
capture nwlsmselect und, clusters(0 2) nolist gname(x)
assert _rc == 198
capture nwlsmselect und, clusters(0 2) nolist gopts(title("x"))
assert _rc == 198
quietly nwlsmselect und, clusters(0 2) seed(2) nsample(200) burnin(200) refit graph
local bst = r(best)
assert e(clusters) == `bst'
graph drop _all

* ---- data in memory unaffected by a plot
quietly {
	clear
	set obs 7
	gen x = _n
}
nwlsmplot, draws(10) name(z1) 
assert _N == 7
summarize x, meanonly
assert r(sum) == 28

* ---- labels, value-labelled legend, noarrows, gray palette, black titles (directed fit of the 18 monks)
quietly nwlsm dir, dim(2) seed(3) nsample(300) burnin(300)
quietly {
	clear
	set obs 18
	gen byte grpn = 1 + (_n > 9)
	label define fac_l 1 "Loyal" 2 "Turks"
	label values grpn fac_l
	label variable grpn "Faction"
	gen byte plain = 1 + mod(_n, 3)
}
local tmpd "$LSM_SCRATCH"
capture program drop lsmt_png
program lsmt_png
	args fn
	quietly checksum "`fn'"
	c_local cks = r(checksum)
end
* numeric label variable, with and without value labels
nwlsmplot, label(grpn) name(lb1)
lsmt_gexists lb1
nwlsmplot, label(plain) name(lb2)
lsmt_gexists lb2
* color(): value labels in the legend (the legend text is in the graph; its file differs from the unlabelled variable's)
nwlsmplot, color(grpn) name(vl1) export("`tmpd'/vl1.png")
lsmt_gexists vl1
quietly label drop fac_l
nwlsmplot, color(grpn) name(vl2) export("`tmpd'/vl2.png")
lsmt_gexists vl2
lsmt_png "`tmpd'/vl1.png"
local c1 `cks'
lsmt_png "`tmpd'/vl2.png"
if "`c1'" == "`cks'" {
	di as error "FAIL: color() legend does not use the value labels"
	exit 9
}
quietly label define fac_l 1 "Loyal" 2 "Turks"
* noarrows works (no arrowheads: a different picture than the default arrows), arrows and noarrows cannot be combined
nwlsmplot, export("`tmpd'/ar1.png") name(ar1)
lsmt_png "`tmpd'/ar1.png"
local c1 `cks'
nwlsmplot, noarrows export("`tmpd'/ar2.png") name(ar2)
lsmt_png "`tmpd'/ar2.png"
if "`c1'" == "`cks'" {
	di as error "FAIL: noarrows has no effect"
	exit 9
}
capture nwlsmplot, arrows noarrows
assert _rc == 198
nwlsmplot, arrowsize(vsmall) name(as1)
lsmt_gexists as1
* colors(): gray and a colour list
nwlsmplot, color(grpn) colors(gray) name(gr1) export("`tmpd'/gr1.png")
lsmt_gexists gr1
nwlsmplot, color(grpn) colors(red blue) name(gr2)
lsmt_gexists gr2
* black titles: the default title is black (the PNG has no navy title pixels: compare with an explicit navy title)
nwlsmplot, noedges name(bt1) export("`tmpd'/bt1.png")
nwlsmplot, noedges title("Latent positions", color(navy)) name(bt2) export("`tmpd'/bt2.png")
lsmt_png "`tmpd'/bt1.png"
local c1 `cks'
lsmt_png "`tmpd'/bt2.png"
if "`c1'" == "`cks'" {
	di as error "FAIL: the default title colour is not distinguishable from navy"
	exit 9
}
graph drop _all
display "PASS test_nwlsm_plots: plot options (arrows, clustermeans, circles, draws, certainty, labels, noaxes, export) and nwlsmselect graph"
