* Goodness of fit (estat gof), posterior predictive check (estat ppcheck) and trace plots (estat mcmcdiag, graph).
*  (1) against the R package latentnet (gof() on three fits, reference files dev/ref/gof_*.csv): the statistics of the OBSERVED network (degree, in- and out-degree,
*      edgewise and dyadwise shared partners, geodesic distances) must be identical cell by cell; the means of the simulated statistics (500
*      networks, one per posterior draw) must agree with latentnet's mean over its 1000 simulated networks within 4 standard errors of the difference
*      (sd of latentnet's simulations, sqrt(1/500 + 1/1000)) in every cell with a positive sd (0.02 in cells where latentnet never saw a count), and correlate with it at least 0.99 across the cells of every statistic;
*  (2) an independent computation in Mata of the statistics of the observed valued network (Poisson, directed) equals the plugin's;
*  (3) the mean total value of the simulated networks equals the sum of the posterior mean expected values e(P) within 3% (all families);
*  (4) results do not depend on threads(), the seed gives the same table, error paths; graphs are drawn.
do cscripts/nwlsm_testlib.do
set linesize 255
global REF "$LSM_REPO/dev/ref"
capture program drop lsmt_gofcase
program lsmt_gofcase
	syntax name(name=case), NET(string) N(integer) DIR(integer) D(integer) GOF(string) [OPTS(string) X(string)]
	lsmt_mknet net_`case', edges("$REF/`net'_edges.csv") n(`n') dir(`dir')
	local opt "`opts'"
	if "`x'" != "" {
		preserve
		import delimited "$REF/`x'", clear varnames(1)
		mkmat x, matrix(__x)
		mata: __x = st_matrix("__x"); st_matrix("__AD", abs(__x * J(1, `n', 1) - J(`n', 1, 1) * __x'))
		restore
		local opt "`opt' edgecov(__AD)"
	}
	quietly nwlsm net_`case', dim(`d') `opt' seed(11) nsample(2000)
	quietly estat gof, nsim(500) seed(5) gof(`gof')
	foreach t of local gof {
		matrix GG_`t' = r(gof_`t')
	}
	local ok 1
	preserve
	import delimited "$REF/gof_`case'.csv", clear varnames(1)
	tempfile ref
	save `ref'
	restore
	foreach t of local gof {
		matrix G = GG_`t'
		preserve
		use `ref', clear
		keep if stat == "`t'"
		local nr = _N
		mkmat obs mean sd, matrix(RR)
		restore
		* the model statistics are given in the order edges, covariates: the same order here
		local nc = rowsof(G)
		if `nc' != `nr' {
			* latentnet's distance has n cells (distance 1 .. n - 1 and Inf); ours the same; model statistics: edges plus covariates
			di as error "`case' `t': `nc' cells here, `nr' in the reference"
			local ok 0
			continue
		}
		local bad 0
		local worst 0
		forvalues i = 1/`nc' {
			if G[`i', 1] != RR[`i', 1] {
				di as error "`case' `t' cell `i': observed " G[`i', 1] " against latentnet " RR[`i', 1]
				local bad 1
			}
			if RR[`i', 3] > 0 {
				local z = abs(G[`i', 3] - RR[`i', 2]) / (RR[`i', 3] * sqrt(1/500 + 1/1000))
				if `z' > `worst' local worst = `z'
				if `z' > 4 {
					di as error "`case' `t' cell `i': simulated mean " G[`i', 3] " against latentnet " RR[`i', 2] " (z = " `z' ")"
					local bad 1
				}
			}
			else if abs(G[`i', 3] - RR[`i', 2]) > 0.02 {
				di as error "`case' `t' cell `i': simulated mean " G[`i', 3] " against latentnet " RR[`i', 2] " (sd 0; rare events: tolerance 0.02)"
				local bad 1
			}
		}
		mata: st_local("cr", strofreal(`nc' > 2 ? correlation((st_matrix("G")[., 3], st_matrix("RR")[., 2]))[1, 2] : 1, "%9.5f"))
		di as text "case `case' `t': observed statistics identical: " cond(`bad', "NO", "yes") "; largest z of the simulated means " %5.2f `worst' "; correlation of the mean vectors " `cr'
		if `bad' | `cr' < 0.99 local ok 0
	}
	if !`ok' {
		di as error "FAIL gof case `case'"
		exit 9
	}
	di as result "PASS gof case `case'"
end
lsmt_gofcase sampson, net(sampson_d2) n(18) dir(1) d(2) gof(idegree odegree espartners dspartners distance)
lsmt_gofcase random40, net(random40_d2) n(40) dir(0) d(2) x(random40_x.csv) gof(degree espartners dspartners distance)
lsmt_gofcase sr40, net(sr40) n(40) dir(1) d(2) opts(sender receiver) x(sr40_x.csv) gof(idegree odegree espartners dspartners distance)
* (2) the observed statistics of a valued network, computed independently in Mata
lsmt_mkvnet net_pois30, y("$REF/pois30_y.csv") dir(1)
quietly nwlsm net_pois30, family(poisson) nsample(300) burnin(500) seed(3)
quietly estat gof, nsim(40) seed(1) gof(idegree odegree espartners dspartners distance model values)
foreach t in idegree odegree espartners dspartners distance model values {
	matrix GV_`t' = r(gof_`t')
}
mata:
mata set matastrict off
void lsmt_obsstats(string scalar Yname, real scalar dir)
{
	real matrix Y, B, Dm, T
	real scalar n, i, j, k, c
	real colvector id, od, hi, ho
	real matrix TT
	real matrix e, dd

	Y = st_matrix(Yname)
	n = rows(Y)
	_diag(Y, 0)
	B = (Y :> 0)
	id = colsum(B)'
	od = rowsum(B)
	hi = J(n, 1, 0)
	ho = J(n, 1, 0)
	for (k = 0; k < n; k++) {
		hi[k + 1] = sum(id :== k)
		ho[k + 1] = sum(od :== k)
	}
	st_matrix("O_ideg", hi)
	st_matrix("O_odeg", ho)
	e = J(n - 1, 1, 0)
	dd = J(n - 1, 1, 0)
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j) continue
		c = sum(B[i, .] :* B[., j]')
		dd[c + 1] = dd[c + 1] + 1
		if (B[i, j]) e[c + 1] = e[c + 1] + 1
	}
	st_matrix("O_esp", e)
	st_matrix("O_dsp", dd)
	Dm = B + (!B) * 1e9
	_diag(Dm, 0)
	for (k = 1; k <= n; k++) {
		TT = Dm[., k] * J(1, n, 1) + J(n, 1, 1) * Dm[k, .]
		Dm = Dm :* (Dm :<= TT) :+ TT :* (Dm :> TT)
	}
	T = J(n, 1, 0)
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j) continue
		if (Dm[i, j] >= 1e8) T[n] = T[n] + 1
		else T[Dm[i, j]] = T[Dm[i, j]] + 1
	}
	st_matrix("O_dist", T)
}
end
preserve
import delimited "$REF/pois30_y.csv", clear varnames(1)
mkmat *, matrix(__Y0)
restore
mata: lsmt_obsstats("__Y0", 1)
foreach pair in idegree:ideg odegree:odeg espartners:esp dspartners:dsp distance:dist {
	gettoken t g : pair, parse(":")
	local g = substr("`g'", 2, .)
	matrix G = GV_`t'
	mata: st_local("dif", strofreal(max(abs(st_matrix("G")[., 1] - st_matrix("O_`g'")))))
	di as text "valued network, observed `t': largest difference to the independent computation " `dif'
	assert `dif' == 0
}
* cutoff(): the structure statistics of a threshold network (tie if the value exceeds the cutoff), Poisson with cutoff(1) and normal with cutoff(0.5), against the
* same independent computation applied to the matrix of the values minus the cutoff
lsmt_mkvnet net_gau30, y("$REF/gau30_y.csv") dir(1)
foreach spec in poisson:pois30:1 normal:gau30:0.5 normal:gau30:-0.7 {
	tokenize "`spec'", parse(":")
	local fam "`1'"
	local nmx "`3'"
	local cut "`5'"
	if "`fam'" == "poisson" {
		quietly nwlsm net_pois30, family(poisson) nsample(300) burnin(500) seed(3)
	}
	else {
		quietly nwlsm net_gau30, family(normal) nsample(200) burnin(300) seed(3)
	}
	quietly estat gof, nsim(30) seed(1) cutoff(`cut') gof(idegree odegree espartners dspartners distance)
	assert abs(r(cutoff) - `cut') < 1e-12
	foreach t in idegree odegree espartners dspartners distance {
		matrix GC_`t' = r(gof_`t')
	}
	preserve
	import delimited "$REF/`nmx'_y.csv", clear varnames(1)
	mkmat *, matrix(__Yc0)
	restore
	mata: st_matrix("__Yc1", st_matrix("__Yc0") :- `cut')
	mata: lsmt_obsstats("__Yc1", 1)
	foreach pair in idegree:ideg odegree:odeg espartners:esp dspartners:dsp distance:dist {
		gettoken t g : pair, parse(":")
		local g = substr("`g'", 2, .)
		matrix G = GC_`t'
		mata: st_local("dif", strofreal(max(abs(st_matrix("G")[., 1] - st_matrix("O_`g'")))))
		di as text "`fam' cutoff(`cut'), observed `t': largest difference to the independent computation " `dif'
		assert `dif' == 0
	}
}
* default statistics with a cutoff after the normal family, cutoff() without a valued family
quietly estat gof, nsim(20) seed(1) cutoff(0)
confirm matrix r(gof_model)
confirm matrix r(gof_values)
capture estat gof, nsim(20) cutoff(abc)
assert _rc == 198
* model statistics: the total value and the covariate sums over the ties
preserve
import delimited "$REF/pois30_x.csv", clear varnames(1)
mkmat x, matrix(__x)
restore
mata: __x = st_matrix("__x"); st_matrix("__AD", abs(__x * J(1, 30, 1) - J(30, 1, 1) * __x'))
quietly nwlsm net_pois30, family(poisson) nsample(200) burnin(300) seed(3) edgecov(__AD)
quietly estat gof, nsim(10) seed(1) gof(model)
matrix G = r(gof_model)
mata: Y = st_matrix("__Y0"); _diag(Y, 0); st_numscalar("__m1", sum(Y)); st_numscalar("__m2", sum(Y :* st_matrix("__AD")))
di as text "model statistics: total value " G[1, 1] " (direct " scalar(__m1) "), sum of the covariate over the ties " G[2, 1] " (direct " scalar(__m2) ")"
assert G[1, 1] == scalar(__m1)
assert abs(G[2, 1] - scalar(__m2)) < 1e-9
* (3) the mean total value of the simulated networks against the sum of e(P)
capture program drop lsmt_totval
program lsmt_totval
	syntax anything(name=case), OPTS(string) NET(string) DIR(integer)
	quietly nwlsm `net', `opts' nsample(600) burnin(1000) seed(4)
	quietly estat gof, nsim(300) seed(2) gof(model)
	matrix G = r(gof_model)
	mata: st_local("sp", strofreal(sum(st_matrix("e(P)")) - sum(diagonal(st_matrix("e(P)"))) , "%20.8f"))
	if !`dir' local sp = `sp' / 2
	di as text "`case': mean total value of the simulated networks " %10.3f G[1, 3] "  sum of the posterior mean expected values " %10.3f `sp'
	assert abs(G[1, 3] - `sp') < 0.03 * abs(`sp') + 1
end
lsmt_mkvnet net_bin30, y("$REF/bin30_y.csv") dir(1)
lsmt_mkvnet net_pois30s, y("$REF/pois30s_y.csv") dir(0)
lsmt_totval poisson, net(net_pois30) opts(family(poisson)) dir(1)
lsmt_totval binomial, net(net_bin30) opts(family(binomial) trials(6) sender receiver) dir(1)
lsmt_totval normal, net(net_gau30) opts(family(normal) dispvar(1) dispdf(3)) dir(1)
lsmt_totval poissonu, net(net_pois30s) opts(family(poisson) sociality) dir(0)
* (4) threads, seed, errors, graphs
quietly nwlsm net_pois30, family(poisson) nsample(300) burnin(500) seed(3)
quietly estat gof, nsim(60) seed(9) threads(1) gof(idegree distance values)
matrix A = r(gof_distance)
matrix AV = r(gof_values)
quietly estat gof, nsim(60) seed(9) threads(4) gof(idegree distance values)
matrix B = r(gof_distance)
matrix BV = r(gof_values)
mata: assert(max(abs(st_matrix("A") - st_matrix("B"))) == 0 & max(abs(st_matrix("AV") - st_matrix("BV"))) == 0)
quietly estat gof, nsim(60) seed(10) gof(idegree distance values)
matrix C = r(gof_distance)
mata: assert(max(abs(st_matrix("A") - st_matrix("C"))) > 0)
capture estat gof, gof(degree)
assert _rc == 198
capture estat gof, gof(foo)
assert _rc == 198
capture estat gof, nsim(1)
assert _rc == 198
capture estat gof, gof(values) threads(-1)
assert _rc == 198
capture estat gof, seed(-5)
assert _rc == 198
quietly nwlsm net_gau30, family(normal) nsample(200) burnin(300) seed(3)
capture estat gof, gof(degree)
assert _rc == 198
capture estat gof, gof(distance)
assert _rc == 198
lsmt_mknet net40, edges("$REF/random40_d2_edges.csv") n(40) dir(0)
quietly nwlsm net40, nsample(300) burnin(500) seed(3)
capture estat gof, gof(idegree)
assert _rc == 198
capture estat gof, gof(values)
assert _rc == 198
capture estat gof, cutoff(0)
assert _rc == 198
estat gof, nsim(40) seed(1) plot name(gofplot1)
graph export "$LSM_SCRATCH/test_gof_plot.png", replace
estat ppcheck, groups(8) graph name(pp1)
assert r(N) == 780
capture estat ppcheck, groups(1)
assert _rc == 198
quietly nwlsm net40, nsample(400) burnin(500) seed(3) chains(3)
estat mcmcdiag, graph name(trace1)
graph export "$LSM_SCRATCH/test_trace.png", replace
quietly nwlsm net_gau30, family(normal) nsample(200) burnin(300) seed(3)
estat mcmcdiag, nolist
matrix D = r(diag)
local rn : rownames D
assert word("`rn'", rowsof(D)) == "dispersion"
ereturn clear
capture estat gof
assert _rc == 301
display "PASS test_nwlsm_gof: estat gof agrees with latentnet's gof(), valued statistics, estat ppcheck and trace plots"
