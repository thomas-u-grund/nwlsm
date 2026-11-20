* The command lines that the dialogs of nwlsmplot, nwlsmselect and estat (nwlsm_estat) build for the graph options run, the graph exists.
* (The assembly of the command string is tested by tools/dlg_assembly_test.py; here the strings are run.)  Run from the repository root.
do cscripts/nwlsm_testlib.do
* (graphics stay on: Stata on Windows cannot export a graph in batch mode with graphics off, r(693))
set linesize 255
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
lsmt_mknet und, edges("$LSM_REPO/dev/ref/clu45_edges.csv") n(45) dir(0)
quietly nwlsm und, dim(2) clusters(2) seed(3) nsample(400) burnin(500)
local out "$LSM_SCRATCH/dlgplots"
capture mkdir "`out'"
capture erase "`out'/g1.png"
nwlsmplot, export(`out'/g1.png) name(d_lp)
lsmt_gexists d_lp
confirm file "`out'/g1.png"
erase "`out'/g1.png"
estat mcmcdiag, graph name(d_gm)
lsmt_gexists d_gm
estat gof, plot name(d_gg) nsim(20) seed(1)
lsmt_gexists d_gg
estat ppcheck, graph name(d_gp)
lsmt_gexists d_gp
* estat influence, graph needs a fit with influence(): its command line is only assembled in the dialog test
nwlsmselect und, clusters(0 2) graph gname(d_gs) gopts(xlabel(none))
lsmt_gexists d_gs
di as text "PASS test_dlg_plots"
