*! version 1.0.0  20nov2026 nwlsmplot
* Copyright (C) 2026 Thomas U. Grund. Part of nwlsm; GNU GPL-3 (see LICENSE).

capture program drop nwlsmplot
program nwlsmplot, rclass
	version 14
	* Stata reads noarrows as the negation of arrows; a real option of its own is needed
	local 0 = subinstr(`"`0'"', "noarrows", "arrowsoff", .)
	syntax [, DIMs(numlist integer min=2 max=2 >=1) MEAN COLor(varname) SIZE(string) LABel(varname) NOEDGES CLUSters ///
		EDGEColor(string) EDGEWidth(string) MSIZE(string) NAME(string) ARRows ARROWSOFF ARROWSize(string) COLORS(string) CLUSTERMeans CIRCles DRAWs(integer 0) CERTainty ///
		LABSize(string) LABPosition(string) NOAXes NOEQual EXPort(string asis) *]
	if "`e(cmd)'" != "nwlsm" error 301
	tempname Z E
	if "`mean'" != "" matrix `Z' = e(Zmean)
	else matrix `Z' = e(Z)
	local n = rowsof(`Z')
	local d = colsof(`Z')
	if "`dims'" == "" {
		if `d' == 1 {
			di as error "the model has one latent dimension; nwlsmplot needs at least two"
			exit 198
		}
		local dims "1 2"
	}
	local d1 : word 1 of `dims'
	local d2 : word 2 of `dims'
	if `d1' > `d' | `d2' > `d' | `d1' == `d2' {
		di as error "dims() must name two different dimensions between 1 and `d'"
		exit 198
	}
	local clsvec
	if "`clusters'" != "" {
		if "`color'" != "" {
			di as error "clusters and color() cannot be combined"
			exit 198
		}
		if "`e(clusters)'" == "" | "`e(clusters)'" == "0" {
			di as error "the model was fitted without latent clusters (option clusters())"
			exit 198
		}
		tempname CLS
		capture matrix `CLS' = e(class)
		if _rc {
			di as error "the classification was not stored (too many actors for a Stata matrix)"
			exit 498
		}
		local clsvec "`CLS'"
		local color "cluster"
	}
	if "`arrows'" != "" & "`arrowsoff'" != "" {
		di as error "arrows and noarrows cannot be combined"
		exit 198
	}
	if "`arrows'" != "" & "`noedges'" != "" {
		di as error "arrows and noedges cannot be combined"
		exit 198
	}
	if `draws' < 0 | `draws' > 1000 {
		di as error "draws() must be between 1 and 1000"
		exit 198
	}
	if "`label'" == "" & (`"`labsize'"' != "" | `"`labposition'"' != "") {
		di as error "labsize() and labposition() need the option label()"
		exit 198
	}
	if ("`clustermeans'" != "" | "`circles'" != "" | "`certainty'" != "") & ("`e(clusters)'" == "" | "`e(clusters)'" == "0") {
		di as error "clustermeans, circles and certainty need a fit with latent clusters (option clusters())"
		exit 198
	}
	if "`certainty'" != "" & "`clusters'" == "" {
		di as error "certainty needs the option clusters (the markers are coloured by cluster, with opacity by the certainty of the membership)"
		exit 198
	}
	_nwsyntax `e(network)', max(1)
	if `nodes' != `n' {
		di as error "the network `e(network)' no longer has the number of actors of the fit"
		exit 198
	}
	local netobj_ `netobj'
	if "`arrows'" != "" & "`directed'" != "true" {
		di as error "arrows: the network `e(network)' is not directed"
		exit 198
	}
	local usearrows = ("`directed'" == "true" & "`arrowsoff'" == "" & "`noedges'" == "")
	* node size: sender, receiver or sociality effects of the fit, or a numeric variable of the data set
	local sizevar
	local sizeeff
	if `"`size'"' != "" {
		if `"`size'"' == "influence" {
			if "`e(influence)'" != "yes" {
				di as error "the fit has no influence model (option influence())"
				exit 198
			}
			local sizeeff influence
		}
		else if inlist(`"`size'"', "sender", "receiver", "sociality") {
			local __ef `"`e(effects)'"'
			if !`: list size in __ef' {
				di as error "the fit has no {bf:`size'} effect"
				exit 198
			}
			local sizeeff `size'
		}
		else {
			capture confirm numeric variable `size'
			if _rc {
				di as error "size(): {bf:`size'} is neither sender, receiver, sociality nor a numeric variable"
				exit 198
			}
			local sizevar `size'
		}
	}
	if ("`color'" != "" & "`clusters'" == "" | "`sizevar'" != "" | "`label'" != "") & _N != `n' {
		di as error "{p 0 0 2}color(), size() and label() take variables with one observation per actor in the order of the actors of {bf:`e(network)'}; the data set has `=_N' observations.{p_end}"
		exit 198
	}
	* label(): a string variable as it is; a numeric variable by its value labels (the numbers where it has none)
	if "`label'" != "" {
		capture confirm string variable `label'
		if _rc {
			tempvar lbs
			local lbvl : value label `label'
			if "`lbvl'" != "" quietly decode `label', generate(`lbs')
			else quietly generate str40 `lbs' = string(`label')
			quietly replace `lbs' = string(`label') if `lbs' == "" & `label' < .
			local labelvar "`label'"
			local label "`lbs'"
		}
	}
	* color(): legend entries by the value labels of the variable, the legend title by its variable label
	local clv
	local cttl "`color'"
	if "`color'" != "" & "`clusters'" == "" {
		quietly levelsof `color', local(clv)
		local ci 0
		foreach v of local clv {
			local ++ci
			local cl_`ci' : label (`color') `v'
			local cl_`ci' = subinstr(`"`cl_`ci''"', `"""', "", .)
		}
		local cvlab : variable label `color'
		if `"`cvlab'"' != "" local cttl = subinstr(`"`cvlab'"', `"""', "", .)
	}
	* collect the data for the plot (before the data in memory are replaced)
	tempname A
	mata: st_matrix("`A'", *(`netobj_'->get_matrix()))
	local sizevec
	if "`sizeeff'" != "" {
		tempname EF
		if "`sizeeff'" == "influence" {
			tempname WW
			matrix `WW' = e(W)
			mata: st_matrix("`EF'", colsum(st_matrix("`WW'"))')
		}
		else matrix `EF' = e(`sizeeff')
		local sizevec "`EF'"
	}
	local G = cond("`e(clusters)'" == "", 0, e(clusters))
	local bil = ("`e(latent)'" == "bilinear")
	tempname QX CMX ZMX
	if "`certainty'" != "" {
		capture matrix `ZMX' = e(Q)
		if _rc {
			di as error "the membership probabilities were not stored (too many actors for a Stata matrix)"
			exit 498
		}
		mata: st_matrix("`QX'", rowmax(st_matrix("`ZMX'")))
	}
	* circles: the radius is the square root of the posterior mean variance of the cluster
	local rad
	if "`clustermeans'" != "" | "`circles'" != "" {
		matrix `CMX' = e(cmeans)
		if "`circles'" != "" {
			tempname VS
			matrix `VS' = e(vsummary)
			forvalues g = 1/`G' {
				local rw = rownumb(`VS', "sigma2_`g'")
				if `rw' >= . {
					di as error "the cluster variances were not found in e(vsummary)"
					exit 498
				}
				local rad `rad' `=sqrt(`VS'[`rw', 1])'
			}
		}
	}
	* the cloud of posterior draws
	local ndr 0
	if `draws' > 0 {
		matrix `ZMX' = e(Zmean)
		matrix `Z' = e(Z)
		mata: st_numscalar("__nwlsm_nd", LsmPlotDraws("`Z'", "`ZMX'", `d1', `d2', `draws', `bil'))
		local ndr = scalar(__nwlsm_nd)
		scalar drop __nwlsm_nd
		if `ndr' < 1 {
			if `ndr' == 0 di as error "draws(): the posterior draws of the fit are not in memory (they stay in memory after nwlsm until the next fit or {bf:clear all})"
			else di as error "draws(): the draws in memory do not belong to the fit in e(); fit the model again"
			exit 498
		}
		if "`mean'" != "" matrix `Z' = e(Zmean)
		else matrix `Z' = e(Z)
	}
	mata: LsmPlotPrep("`A'", "`Z'", `d1', `d2', `=("`directed'" == "true")', "`=cond("`clusters'" != "", "", "`color'")'", "`clsvec'", "`sizevar'", "`sizevec'", "`label'", `=("`noedges'" == "")', `usearrows')
	mata: st_local("m", strofreal(rows(__nwlsm_ed)))
	local ncirc = cond("`circles'" != "", 101 * `G', 0)
	local nmean = cond("`clustermeans'" != "", `G', 0)
	local nd = `ndr' * `n'
	preserve
	quietly {
		clear
		set obs `=`n' + `m' + `nd' + `ncirc' + `nmean''
		gen double x = .
		gen double y = .
		gen double cg = .
		gen double w = .
		gen str40 lb = ""
		gen double ex1 = .
		gen double ey1 = .
		gen double ex2 = .
		gen double ey2 = .
		gen double qb = .
		gen double dx = .
		gen double dy = .
		gen double da = .
		gen double dg = .
		gen double cx = .
		gen double cy = .
		gen double cgr = .
		gen double mx = .
		gen double my = .
		gen double mg = .
		mata: LsmPlotStore()
		if "`certainty'" != "" mata: st_store((1::`n'), "qb", 1 :+ floor(3.999 :* rowmax((st_matrix("`QX'") :- 1/`G', J(`n', 1, 0))) :/ (1 - 1/`G')))
		if `nd' > 0 {
			mata: LsmPlotDrawStore(`n' + `m' + 1)
			replace dg = cg[da] in `=`n' + `m' + 1'/`=`n' + `m' + `nd''
		}
		local o = `n' + `m' + `nd'
		if `ncirc' > 0 {
			local k 0
			forvalues g = 1/`G' {
				local r : word `g' of `rad'
				forvalues t = 0/100 {
					local ++k
					replace cx = `CMX'[`g', `d1'] + `r' * cos(2 * c(pi) * `t' / 100) in `=`o' + `k''
					replace cy = `CMX'[`g', `d2'] + `r' * sin(2 * c(pi) * `t' / 100) in `=`o' + `k''
					replace cgr = `g' in `=`o' + `k''
				}
			}
		}
		local o = `o' + `ncirc'
		if `nmean' > 0 {
			forvalues g = 1/`G' {
				replace mx = `CMX'[`g', `d1'] in `=`o' + `g''
				replace my = `CMX'[`g', `d2'] in `=`o' + `g''
				replace mg = `g' in `=`o' + `g''
			}
		}
	}
	local pal navy maroon forest_green dkorange teal purple gold cyan olive magenta brown gs8 ltblue red green
	* colors(): a list of colours for the groups; gray (or grey): groups differ by gray value and by marker (filled and hollow)
	local syms
	local gray 0
	if `"`colors'"' != "" {
		if inlist(lower(`"`colors'"'), "gray", "grey") {
			local pal gs2 gs10 gs5 gs12 gs7 gs3 gs9 gs13 gs4 gs11 gs6 gs1 gs8 gs14 gs2
			local syms O Dh T Sh X D Th S Oh + Dh Sh T O D
			local gray 1
		}
		else local pal `"`colors'"'
	}
	local plots
	local nplot 0
	if "`noedges'" == "" & `m' > 0 {
		* dense networks: lighter lines and smaller arrowheads
		local dense = (`m' > 400)
		if `"`edgecolor'"' == "" local edgecolor = cond(`dense', "gs12%35", "gs11%60")
		if `"`edgewidth'"' == "" local edgewidth "vthin"
		if `"`arrowsize'"' == "" local arrowsize = cond(`dense', "vsmall", "small")
		if `usearrows' local plots `"(pcarrow ey1 ex1 ey2 ex2, lcolor(`edgecolor') mcolor(`edgecolor') lwidth(`edgewidth') barbsize(`arrowsize'))"'
		else local plots `"(pcspike ey1 ex1 ey2 ex2, lcolor(`edgecolor') lwidth(`edgewidth'))"'
		local ++nplot
	}
	* levels of the colour variable
	local lv
	if "`color'" != "" {
		quietly levelsof cg, local(lv)
		local nl : word count `lv'
		if `nl' > 15 {
			di as error "color(): the variable has more than 15 distinct values; recode it into groups"
			exit 198
		}
	}
	* cloud of posterior draws (under the actors)
	if `ndr' > 0 {
		if "`color'" != "" {
			local g 0
			foreach v of local lv {
				local ++g
				local pi = cond("`clusters'" != "", `v', `g')
				local plots `"`plots' (scatter dy dx if dg == `v', msymbol(point) msize(small) mcolor(`: word `pi' of `pal''%12))"'
				local ++nplot
			}
		}
		else {
			local plots `"`plots' (scatter dy dx, msymbol(point) msize(small) mcolor(navy%12))"'
			local ++nplot
		}
	}
	* variance circles and cluster means
	if `ncirc' > 0 {
		forvalues g = 1/`G' {
			local cc = cond("`clusters'" != "", "`: word `g' of `pal''", "black")
			local plots `"`plots' (line cy cx if cgr == `g', lcolor(`cc') lpattern(dash) lwidth(thin))"'
			local ++nplot
		}
	}
	if `nmean' > 0 {
		forvalues g = 1/`G' {
			local cc = cond("`clusters'" != "", "`: word `g' of `pal''", "black")
			local plots `"`plots' (scatter my mx if mg == `g', msymbol(plus) msize(large) mcolor(`cc'))"'
			local ++nplot
		}
	}
	local mopt
	if `"`msize'"' == "" local msize "medium"
	if `"`labsize'"' == "" local labsize "vsmall"
	if `"`labposition'"' == "" local labposition "auto"
	local autolab 0
	if `"`label'"' != "" {
		if `"`labposition'"' == "auto" {
			if `n' <= 150 local autolab 1
			else local labposition "1"
		}
		if `autolab' local mopt `"mlabel(lb) mlabsize(`labsize') mlabvposition(lp)"'
		else local mopt `"mlabel(lb) mlabsize(`labsize') mlabposition(`labposition')"'
	}
	local wt
	if "`sizevar'`sizeeff'" != "" local wt "[aweight = w]"
	else local mopt `"`mopt' msize(`msize')"'
	local legend
	if "`color'" != "" {
		local pal2 `"`pal'"'
		local g 0
		local order
		foreach v of local lv {
			local ++g
			local pi = cond("`clusters'" != "", `v', `g')
			local col : word `pi' of `pal2'
			local sym
			local lcl `col'
			if `gray' {
				local sym "msymbol(`: word `pi' of `syms'')"
				local lcl black
			}
			local lab "`v'"
			if "`clusters'" == "" {
				local ci : list posof "`v'" in clv
				if `ci' > 0 local lab `"`cl_`ci''"'
			}
			if "`certainty'" != "" {
				local last
				local bi 0
				foreach op in 30 50 75 100 {
					local ++bi
					quietly count if cg == `v' & qb == `bi'
					if r(N) > 0 {
						local ++nplot
						local last `nplot'
						local plots `"`plots' (scatter y x `wt' if cg == `v' & qb == `bi', mcolor(`col'%`op') mlabcolor(`lcl') `sym' `mopt')"'
					}
				}
				local order `"`order' `last' "`lab'""'
			}
			else {
				local ++nplot
				local plots `"`plots' (scatter y x `wt' if cg == `v', mcolor(`col'`=cond(`gray', "", "%75")') mlabcolor(`lcl') `sym' `mopt')"'
				local order `"`order' `nplot' "`lab'""'
			}
		}
		local legend `"legend(order(`order') position(3) cols(1) title(`"`cttl'"', size(small) color(black)))"'
	}
	else {
		local plots `"`plots' (scatter y x `wt', mcolor(navy%75) mlabcolor(navy) `mopt')"'
		local legend "legend(off)"
	}
	local ttl `"title("Latent positions", color(black))"'
	if strpos(`"`options'"', "title(") local ttl
	local axopt `"xtitle("dimension `d1'", color(black)) ytitle("dimension `d2'", color(black)) xlabel(, grid) ylabel(, grid)"'
	if "`noaxes'" != "" local axopt `"xtitle("") ytitle("") xscale(off) yscale(off) xlabel(none) ylabel(none)"'
	* equal scale: the same data range (plus the same padding) on both axes' range lengths, aspect ratio = y range / x range, so that one unit
	* has the same length on x and y (distances are the content of a latent space; circles stay round)
	local aspopt "aspectratio(1)"
	if "`noequal'" == "" {
		quietly summarize x, meanonly
		local xlo = r(min)
		local xhi = r(max)
		quietly summarize y, meanonly
		local ylo = r(min)
		local yhi = r(max)
		foreach v in ex1 ex2 dx cx mx {
			quietly summarize `v', meanonly
			if r(N) > 0 {
				local xlo = min(`xlo', r(min))
				local xhi = max(`xhi', r(max))
			}
		}
		foreach v in ey1 ey2 dy cy my {
			quietly summarize `v', meanonly
			if r(N) > 0 {
				local ylo = min(`ylo', r(min))
				local yhi = max(`yhi', r(max))
			}
		}
		local xr = max(`xhi' - `xlo', 1e-6)
		local yr = max(`yhi' - `ylo', 1e-6)
		local pad = 0.06 * max(`xr', `yr')
		local xlo = `xlo' - `pad'
		local xhi = `xhi' + `pad'
		local ylo = `ylo' - `pad'
		local yhi = `yhi' + `pad'
		local asp = (`yhi' - `ylo') / (`xhi' - `xlo')
		* tick labels inside the range, with one step on both axes (a label outside the range would widen the axis and break the equal scale)
		local raw = max(`xhi' - `xlo', `yhi' - `ylo') / 6
		local mag = 10^floor(log10(`raw'))
		local fr = `raw' / `mag'
		local step = `mag' * cond(`fr' <= 1.5, 1, cond(`fr' <= 3.5, 2, cond(`fr' <= 7.5, 5, 10)))
		local xl
		forvalues k = `=ceil(`xlo' / `step')'/`=floor(`xhi' / `step')' {
			local xl `xl' `=`k' * `step''
		}
		local yl
		forvalues k = `=ceil(`ylo' / `step')'/`=floor(`yhi' / `step')' {
			local yl `yl' `=`k' * `step''
		}
		local aspopt "aspectratio(`asp') xscale(range(`xlo' `xhi')) yscale(range(`ylo' `yhi'))"
		if "`xl'" != "" & "`yl'" != "" & "`noaxes'" == "" local axopt `"xtitle("dimension `d1'", color(black)) ytitle("dimension `d2'", color(black)) xlabel(`xl', grid) ylabel(`yl', grid)"'
	}
	if strpos(`"`options'"', "aspectratio(") | strpos(`"`options'"', "aspect(") local aspopt
	* automatic label positions: for every actor the clock position (around the marker) with the fewest other markers and already placed labels under its label
	if `autolab' {
		quietly {
			summarize x in 1/`n', meanonly
			local axlo = r(min)
			local axhi = r(max)
			summarize y in 1/`n', meanonly
			local aylo = r(min)
			local ayhi = r(max)
			if "`noequal'" == "" {
				local axlo = min(`axlo', `xlo')
				local axhi = max(`axhi', `xhi')
				local aylo = min(`aylo', `ylo')
				local ayhi = max(`ayhi', `yhi')
			}
			else {
				local pd = 0.06 * max(`axhi' - `axlo', `ayhi' - `aylo', 1e-6)
				local axlo = `axlo' - `pd'
				local axhi = `axhi' + `pd'
				local aylo = `aylo' - `pd'
				local ayhi = `ayhi' + `pd'
			}
			local lfac = cond("`labsize'" == "vsmall", 1, cond("`labsize'" == "small", 1.3, cond("`labsize'" == "medsmall", 1.6, cond("`labsize'" == "medium", 2, 1.2))))
			generate double __nx = (x - `axlo') / (`axhi' - `axlo') in 1/`n'
			generate double __ny = (y - `aylo') / (`ayhi' - `aylo') in 1/`n'
			generate double __lw = (0.0125 * max(1, strlen(lb)) * `lfac') in 1/`n'
			local lh = 0.042 * `lfac'
			generate double __bx = . 
			generate double __by = . 
			generate double __bw = . 
			generate byte lp = 3 in 1/`n'
			local cpos "3 9 12 6 1 11 5 7"
			local csx "1 -1 0 0 1 -1 1 -1"
			local csy "0 0 1 -1 1 1 -1 -1"
			forvalues i = 1/`n' {
				local xi = __nx[`i']
				local yi = __ny[`i']
				local wi = __lw[`i']
				local best = 1e9
				forvalues c = 1/8 {
					local sx : word `c' of `csx'
					local sy : word `c' of `csy'
					local bxc = `xi' + `sx' * (`wi' / 2 + 0.008)
					local byc = `yi' + `sy' * (`lh' / 2 + 0.008)
					count if abs(__nx - `bxc') < `wi' / 2 + 0.004 & abs(__ny - `byc') < `lh' / 2 + 0.004 & _n != `i' & _n <= `n'
					local pen = 2 * r(N) + 0.01 * `c'
					count if abs(__bx - `bxc') < (`wi' + __bw) / 2 & abs(__by - `byc') < `lh'
					local pen = `pen' + 3 * r(N)
					if `bxc' - `wi' / 2 < 0 | `bxc' + `wi' / 2 > 1 | `byc' - `lh' / 2 < 0 | `byc' + `lh' / 2 > 1 local pen = `pen' + 1
					if `pen' < `best' {
						local best = `pen'
						local bp : word `c' of `cpos'
						local bbx = `bxc'
						local bby = `byc'
					}
				}
				replace lp = `bp' in `i'
				replace __bx = `bbx' in `i'
				replace __by = `bby' in `i'
				replace __bw = `wi' in `i'
			}
		}
	}
	twoway `plots', `axopt' `ttl' `aspopt' `legend' ///
		`=cond("`name'" != "", "name(`name', replace)", "")' `options'
	if `"`export'"' != "" {
		gettoken __fn __ex : export, parse(",")
		local __ex = trim(substr(`"`__ex'"', 2, .))
		graph export `__fn', replace `__ex'
	}
	if "`noequal'" == "" {
		return scalar aspect = `asp'
		return scalar xmin = `xlo'
		return scalar xmax = `xhi'
		return scalar ymin = `ylo'
		return scalar ymax = `yhi'
	}
	restore
end
