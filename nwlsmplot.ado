*! version 0.3.0  19nov2026 nwlsmplot
* Copyright (C) 2026 Thomas U. Grund. Part of nwlsm; GNU GPL-3 (see LICENSE).

capture program drop nwlsmplot
program nwlsmplot
	version 14
	syntax [, DIMs(numlist integer min=2 max=2 >=1) MEAN COLor(varname) SIZE(string) LABel(varname string) NOEDGES ///
		EDGEColor(string) EDGEWidth(string) MSIZE(string) NAME(string) *]
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
	_nwsyntax `e(network)', max(1)
	if `nodes' != `n' {
		di as error "the network `e(network)' no longer has the number of actors of the fit"
		exit 198
	}
	local netobj_ `netobj'
	* node size: sender, receiver or sociality effects of the fit, or a numeric variable of the data set
	local sizevar
	local sizeeff
	if `"`size'"' != "" {
		if inlist(`"`size'"', "sender", "receiver", "sociality") {
			if !`: list size in local(e(effects))' & `"`size'"' != "" {
				local __ef `"`e(effects)'"'
				if !`: list size in __ef' {
					di as error "the fit has no {bf:`size'} effect"
					exit 198
				}
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
	if ("`color'" != "" | "`sizevar'" != "" | "`label'" != "") & _N != `n' {
		di as error "{p 0 0 2}color(), size() and label() take variables with one observation per actor in the order of the actors of {bf:`e(network)'}; the data set has `=_N' observations.{p_end}"
		exit 198
	}
	* collect the data for the plot (before the data in memory are replaced)
	tempname A
	mata: st_matrix("`A'", *(`netobj_'->get_matrix()))
	local sizevec
	if "`sizeeff'" != "" {
		tempname EF
		matrix `EF' = e(`sizeeff')
		local sizevec "`EF'"
	}
	mata: LsmPlotPrep("`A'", "`Z'", `d1', `d2', `=("`directed'" == "true")', "`color'", "`sizevar'", "`sizevec'", "`label'", `=("`noedges'" == "")')
	local m = rows(st_matrix("`A'")) * 0
	mata: st_local("m", strofreal(rows(__nwlsm_ed)))
	preserve
	quietly {
		clear
		set obs `=`n' + `m''
		gen double x = .
		gen double y = .
		gen double cg = .
		gen double w = .
		gen str40 lb = ""
		gen double ex1 = .
		gen double ey1 = .
		gen double ex2 = .
		gen double ey2 = .
		mata: LsmPlotStore()
	}
	local plots
	if "`noedges'" == "" & `m' > 0 {
		if `"`edgecolor'"' == "" local edgecolor "gs11%60"
		if `"`edgewidth'"' == "" local edgewidth "vthin"
		local plots `"(pcspike ey1 ex1 ey2 ex2, lcolor(`edgecolor') lwidth(`edgewidth'))"'
	}
	local mopt
	if `"`msize'"' == "" local msize "medium"
	if `"`label'"' != "" local mopt `"mlabel(lb) mlabsize(vsmall) mlabposition(1)"'
	local wt
	if "`sizevar'`sizeeff'" != "" local wt "[aweight = w]"
	else local mopt `"`mopt' msize(`msize')"'
	local legend
	if "`color'" != "" {
		quietly levelsof cg, local(lv)
		local nl : word count `lv'
		if `nl' > 15 {
			di as error "color(): the variable has more than 15 distinct values; recode it into groups"
			exit 198
		}
		local g 0
		local order
		local nplot = cond("`noedges'" == "" & `m' > 0, 1, 0)
		foreach v of local lv {
			local ++g
			local ++nplot
			local lab "`v'"
			local plots `"`plots' (scatter y x `wt' if cg == `v', mcolor(p`g'%80) `mopt')"'
			local order `"`order' `nplot' "`lab'""'
		}
		local legend `"legend(order(`order') position(3) cols(1) title("`color'", size(small)))"'
	}
	else {
		local plots `"`plots' (scatter y x `wt', mcolor(navy%75) `mopt')"'
		local legend "legend(off)"
	}
	local ttl "Latent positions"
	twoway `plots', xtitle("dimension `d1'") ytitle("dimension `d2'") title("`ttl'") aspectratio(1) xlabel(, grid) ylabel(, grid) `legend' ///
		`=cond("`name'" != "", "name(`name', replace)", "")' `options'
	restore
end
