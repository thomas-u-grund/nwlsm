* BIC of the latent position model (estat bic, nwlsmselect) against the R package latentnet (bic.ergmm; reference results in dev/ref, made by dev/ref_latentnet_bic.R).
* Part 1 (exact): the BIC components computed from latentnet's own minimum Kullback-Leibler positions and effects (Y part: conditional maximum likelihood
*   of the coefficients, Z part: normal or EM mixture fit, effects parts) equal latentnet's numbers (Sampson monks, directed; a planted-cluster network, undirected).
* Part 2 (statistical): nwlsm fits of the same networks give components close to latentnet's, the same ranking over the number of clusters.
do cscripts/nwlsm_testlib.do
global LSM_OK 1
capture program drop lsmt_bicexact
program lsmt_bicexact
	syntax name(name=case), EDGES(string) N(integer) DIR(integer) G(integer) [EFF(integer 0)]
	local ref "$LSM_REPO/dev/ref"
	preserve
	import delimited "`ref'/`edges'_edges.csv", clear varnames(1)
	mata: lsmt_edges(`n', `dir')
	import delimited "`ref'/`case'_mklZ.csv", clear varnames(1)
	mkmat _all, matrix(ZZ)
	import delimited "`ref'/`case'_bic.csv", clear varnames(1)
	mkmat y z sender receiver overall p edges, matrix(RB)
	local act "0 0 0"
	capture confirm file "`ref'/`case'_mkleff.csv"
	if !_rc {
		import delimited "`ref'/`case'_mkleff.csv", clear varnames(1)
		mkmat sender receiver, matrix(EF)
		local act "1 1 0"
	}
	restore
	mata: st_local("ok", strofreal(lsmt_bicx(`n', `dir', `g', "`act'")))
	if `ok' == 0 {
		di as error "FAIL `case'"
		global LSM_OK 0
	}
	else di as result "PASS `case'"
end
mata:
mata set matastrict off
real scalar lsmt_bicx(real scalar n, real scalar dir, real scalar G, string scalar acts)
{
	real matrix y, Z, E, RB
	real rowvector r, act
	real scalar ok, sbf, rbf

	y = st_matrix("__lsmt_A")
	Z = st_matrix("ZZ")
	RB = st_matrix("RB")
	act = strtoreal(tokens(acts))
	sbf = RB[1, 3]
	rbf = RB[1, 4]
	E = J(n, 3, 0)
	if (act[1]) E[., 1..2] = st_matrix("EF")
	r = LsmBic(y, dir, J(0, n, .), Z, E, act, G, RB[1, 7])
	printf("  nwlsm: Y %10.4f  Z %10.4f  sender %8.4f  receiver %8.4f  overall %10.4f\n", r[1], r[2], r[3], r[4], r[6])
	printf("  R:     Y %10.4f  Z %10.4f  sender %8.4f  receiver %8.4f  overall %10.4f\n", RB[1, 1], RB[1, 2], sbf, rbf, RB[1, 5])
	// the Y and effects parts are exact; the mixture part equals latentnet's unless the EM runs end in different local maxima (ours then never worse)
	ok = (abs(r[1] - RB[1, 1]) < 1e-3) & (abs(r[3] - sbf) < 1e-3) & (abs(r[4] - rbf) < 1e-3)
	if (G == 0) ok = ok & (abs(r[2] - RB[1, 2]) < 1e-3)
	else ok = ok & (r[2] < RB[1, 2] + 1e-3) & (r[2] > RB[1, 2] - 0.05 * abs(RB[1, 2]))
	return(ok)
}
end
lsmt_bicexact bicsamp_g0, edges(sampson_d2) n(18) dir(1) g(0)
lsmt_bicexact bicsamp_g2, edges(sampson_d2) n(18) dir(1) g(2)
lsmt_bicexact bicsamp_g3, edges(sampson_d2) n(18) dir(1) g(3)
lsmt_bicexact bicsamp_sr, edges(sampson_d2) n(18) dir(1) g(0)
lsmt_bicexact bicclu45_g0, edges(clu45) n(45) dir(0) g(0)
lsmt_bicexact bicclu45_g2, edges(clu45) n(45) dir(0) g(2)
lsmt_bicexact bicclu45_g3, edges(clu45) n(45) dir(0) g(3)
lsmt_bicexact bicclu45_g4, edges(clu45) n(45) dir(0) g(4)
* Part 2: complete fits against latentnet's overall BIC (Monte Carlo and MKL differences: within 1.5 per cent), ranking over the number of clusters
capture program drop lsmt_bicfit
program lsmt_bicfit
	syntax name(name=case), EDGES(string) N(integer) DIR(integer) G(string) SEEDS(integer)
	local ref "$LSM_REPO/dev/ref"
	lsmt_mknet net_`case', edges("`ref'/`edges'_edges.csv") n(`n') dir(`dir')
	local ok 1
	foreach g of numlist `g' {
		preserve
		import delimited "`ref'/`case'_g`g'_bic.csv", clear varnames(1)
		local rb = overall[1]
		restore
		if `g' == 0 quietly nwlsm net_`case', dim(2) seed(`seeds')
		else quietly nwlsm net_`case', dim(2) clusters(`g') seed(`seeds')
		quietly estat bic
		local d = abs(r(bic) - `rb') / `rb'
		di as text "case `case' G=`g': nwlsm BIC " %9.3f r(bic) " latentnet " %9.3f `rb' "  relative difference " %6.4f `d'
		if `d' > 0.015 local ok 0
		if abs(r(bic) - (r(bic_y) + r(bic_z))) > 1e-8 local ok 0
	}
	nwlsmselect net_`case', clusters(0 2 3) dim(2) seed(`seeds')
	matrix TT = r(table)
	if r(best) != 3 local ok 0
	if !`ok' {
		di as error "FAIL fit `case'"
		global LSM_OK 0
	}
	else di as result "PASS fit `case'"
end
lsmt_bicfit bicsamp, edges(sampson_d2) n(18) dir(1) g(0 2 3) seeds(1)
lsmt_bicfit bicclu45, edges(clu45) n(45) dir(0) g(0 2 3) seeds(5)
* sender and receiver effects: their MKL estimates are weakly identified; both packages can end in an inferior local optimum (effects of about 0.65 instead of
* 0.27 standard deviation; latentnet: 3 of 12 seeds, BIC 455 to 602 against 398 to 402). With several MKL starts the best optimum is the usual one and the BIC equals latentnet's typical value.
foreach s in 2 3 {
	quietly nwlsm net_bicsamp, dim(2) sender receiver seed(`s') mklstarts(8)
	quietly estat bic, nolist
	di as text "sender and receiver, mklstarts(8), seed `s': BIC " %8.2f r(bic) " (latentnet 399.2; typical 397.9 to 402)"
	if abs(r(bic) - 399.235) > 6 | r(bic_sender) > 12 {
		di as error "FAIL sender and receiver BIC"
		global LSM_OK 0
	}
}
* options of estat bic, results, errors
quietly nwlsm net_bicclu45, dim(2) clusters(3) seed(5)
estat bic, effobs(dyads) nolist
assert abs(r(effobs) - 45 * 44 / 2) < 1e-9
estat bic, effobs(actors) nolist
assert r(effobs) == 45
assert r(clusters_used) == 3
estat bic, nolist
assert abs(r(effobs) - 182) < 1e-9
capture estat bic, effobs(zzz)
assert _rc == 198
quietly nwlsm net_bicclu45, dim(2) bilinear seed(5)
capture estat bic
assert _rc == 198
capture nwlsmselect net_bicclu45, clusters(1 2)
assert _rc == 198
capture nwlsmselect net_bicclu45, clusters(2)
assert _rc == 198
quietly nwlsm net_bicclu45, dim(2) family(poisson) seed(5)
capture estat bic
assert _rc == 198
if $LSM_OK display "PASS test_nwlsm_bic: BIC components equal latentnet's, fits close, ranking over the clusters"
else exit 9
