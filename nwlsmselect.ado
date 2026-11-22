*! version 1.0.6  22nov2026 nwlsmselect
* Copyright (C) 2026 Thomas U. Grund. Part of nwlsm; GNU GPL-3 (see LICENSE).

capture program drop nwlsmselect
program nwlsmselect, rclass
	version 14
	syntax anything(name=netname) [, CLUsters(numlist integer >=0 <=6) EFFobs(string) REFit NOList GRAPH GNAME(string) GOPts(string asis) *]
	_nwlsm_libcheck
	if "`clusters'" == "" local clusters "0 2 3 4"
	local clusters : list uniq clusters
	if `: word count `clusters'' < 2 {
		di as error "clusters() needs at least two values (0 stands for no clusters)"
		exit 198
	}
	foreach g of local clusters {
		if `g' == 1 {
			di as error "clusters() may contain 0 (no clusters) and 2 to 6, not 1"
			exit 198
		}
	}
	if strpos(`"`options'"', "influence(") | strpos(`"`options'"', "bilinear") | strpos(`"`options'"', "family(") {
		di as error "options influence(), bilinear and family() cannot be used with nwlsmselect (binary networks, Euclidean distance)"
		exit 198
	}
	local nm : word count `clusters'
	tempname T
	matrix `T' = J(`nm', 8, .)
	local i 0
	foreach g of local clusters {
		local ++i
		if `g' == 0 quietly nwlsm `netname', `options'
		else quietly nwlsm `netname', `options' clusters(`g')
		if "`effobs'" == "" quietly estat bic, nolist
		else quietly estat bic, effobs(`effobs') nolist
		matrix `T'[`i', 1] = `g'
		matrix `T'[`i', 2] = r(bic)
		matrix `T'[`i', 3] = r(bic_y)
		matrix `T'[`i', 4] = r(bic_z)
		matrix `T'[`i', 5] = r(bic_sender) + r(bic_receiver) + r(bic_sociality)
		matrix `T'[`i', 6] = e(ll_mean)
		matrix `T'[`i', 7] = r(clusters_used)
	}
	* differences to the smallest BIC, the order by BIC
	local bi 1
	forvalues i = 2/`nm' {
		if `T'[`i', 2] < `T'[`bi', 2] local bi `i'
	}
	forvalues i = 1/`nm' {
		matrix `T'[`i', 8] = `T'[`i', 2] - `T'[`bi', 2]
	}
	matrix colnames `T' = clusters BIC BIC_Y BIC_Z BIC_effects loglik used delta
	local bestG = `T'[`bi', 1]
	if "`nolist'" == "" {
		di
		di as text "{p 0 2 2}Selection of the number of latent clusters by BIC (as latentnet's bic.ergmm), network {bf:`netname'}{p_end}"
		di as text "{hline 72}"
		di as text %10s "Clusters" %12s "BIC" %12s "BIC (Y)" %12s "BIC (Z)" %12s "Mean loglik" %12s "BIC - min"
		di as text "{hline 72}"
		local done " "
		forvalues k = 1/`nm' {
			local r 0
			forvalues i = 1/`nm' {
				if strpos("`done'", " `i' ") == 0 {
					if `r' == 0 local r `i'
					else if `T'[`i', 2] < `T'[`r', 2] local r `i'
				}
			}
			local done "`done'`r' "
			di as result %10.0f `T'[`r', 1] %12.3f `T'[`r', 2] %12.3f `T'[`r', 3] %12.3f `T'[`r', 4] %12.3f `T'[`r', 6] %12.3f `T'[`r', 8] as text cond(`r' == `bi', "  <- best", "")
		}
		di as text "{hline 72}"
		di as text "{p 0 0 2}Sorted by BIC, smaller is better; 0 clusters is the model without clusters. The BIC compares fits to the same network with the same options.{p_end}"
	}
	if "`gname'" != "" | `"`gopts'"' != "" {
		if "`graph'" == "" {
			di as error "gname() and gopts() need the option graph"
			exit 198
		}
	}
	if "`graph'" != "" {
		* BIC and its parts against the number of clusters (0 = no clusters)
		preserve
		quietly {
			clear
			svmat double `T', names(col)
			sort clusters
			summarize BIC_effects, meanonly
			local haseff = (r(max) != 0 | r(min) != 0)
			local xl
			forvalues i = 1/`nm' {
				local v = clusters[`i']
				local xl `"`xl' `v'"'
			}
			if clusters[1] == 0 local xl : subinstr local xl " 0" `" 0 "none""', count(local junk)
		}
		local pl `"(connected BIC clusters, lcolor(black) mcolor(black) msymbol(O) lwidth(medthick)) (connected BIC_Y clusters, lcolor(navy) mcolor(navy) msymbol(D) lpattern(dash)) (connected BIC_Z clusters, lcolor(maroon) mcolor(maroon) msymbol(T) lpattern(dash))"'
		local lg `"1 "BIC" 2 "BIC of the ties (Y)" 3 "BIC of the positions (Z)""'
		if `haseff' {
			local pl `"`pl' (connected BIC_effects clusters, lcolor(forest_green) mcolor(forest_green) msymbol(S) lpattern(dash))"'
			local lg `"`lg' 4 "BIC of the random effects""'
		}
		local gtt `"title("BIC by number of latent clusters", color(black))"'
		if strpos(`"`gopts'"', "title(") local gtt
		twoway `pl', xlabel(`xl') xtitle("Number of clusters (none = no clusters)", color(black)) ytitle("BIC", color(black)) xline(`bestG', lpattern(shortdash) lcolor(gs10)) ///
			legend(order(`lg') position(6) rows(2) size(small)) `gtt' note("Dashed vertical line: smallest BIC (`bestG' clusters). Smaller is better.", size(vsmall)) ///
			`=cond("`gname'" != "", "name(`gname', replace)", "")' `gopts'
		restore
	}
	if "`refit'" != "" {
		if `bestG' == 0 quietly nwlsm `netname', `options'
		else quietly nwlsm `netname', `options' clusters(`bestG')
	}
	return matrix table = `T'
	return scalar best = `bestG'
end

