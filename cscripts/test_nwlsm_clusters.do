* Latent clusters (clusters(G)) against the planted clusters and the R package latentnet (reference results in dev/ref, made by dev/ref_latentnet.R):
* an undirected network with three well separated planted clusters (n = 45) and a directed network with two (n = 36, one covariate).
* Criteria (four seeds per case, default sampler settings):
*  - adjusted Rand index of the posterior classification (most probable cluster) with the planted clusters at least 0.9 in every run, and with latentnet's classification at least 0.9;
*  - posterior mean of every coefficient (average over the seeds) within 3.5 Monte Carlo standard errors of the latentnet value, the standard error combining the
*    batch-means error of our average and the difference between two latentnet chains;
*  - posterior membership probabilities: rows sum to 1, the classification equals the row maximum;
*  - the same classification quality with chains(3) and threads(2).
do cscripts/nwlsm_testlib.do
global NS 4
mata:
mata set matastrict off
real scalar lsmt_ari(real colvector a, real colvector b)
{
	real matrix T
	real scalar n, i, j, sij, sa, sb, ex, mx
	real colvector ua, ub

	ua = uniqrows(a)
	ub = uniqrows(b)
	T = J(rows(ua), rows(ub), 0)
	for (i = 1; i <= rows(a); i++) T[selectindex(ua :== a[i]), selectindex(ub :== b[i])] = T[selectindex(ua :== a[i]), selectindex(ub :== b[i])] + 1
	n = rows(a)
	sij = sum(T :* (T :- 1) / 2)
	sa = sum(rowsum(T) :* (rowsum(T) :- 1) / 2)
	sb = sum(colsum(T)' :* (colsum(T)' :- 1) / 2)
	ex = sa * sb / (n * (n - 1) / 2)
	mx = (sa + sb) / 2
	if (mx == ex) return(1)
	return((sij - ex) / (mx - ex))
}
end
mata:
real scalar lsmt_cls_bad(real matrix Q, real colvector CL)
{
	real scalar i, b

	b = 0
	for (i = 1; i <= rows(Q); i++) if (Q[i, CL[i]] < max(Q[i, .]) - 1e-12) b++
	return(b)
}
end
capture program drop lsmt_clu
program lsmt_clu
	syntax name(name=case), N(integer) DIR(integer) G(integer) [X]
	local ref "$LSM_REPO/dev/ref"
	lsmt_mknet net_`case', edges("`ref'/`case'_edges.csv") n(`n') dir(`dir')
	preserve
	import delimited "`ref'/`case'_truth.csv", clear varnames(1)
	mkmat truth, matrix(TR)
	import delimited "`ref'/`case'_class.csv", clear varnames(1)
	mkmat class, matrix(LC)
	import delimited "`ref'/`case'_beta.csv", clear varnames(1)
	mkmat mean sd, matrix(RB)
	import delimited "`ref'/`case'_c2_beta.csv", clear varnames(1)
	mkmat mean, matrix(RB2)
	restore
	local opt
	if "`x'" != "" {
		preserve
		import delimited "`ref'/`case'_x.csv", clear varnames(1)
		mkmat x, matrix(__x)
		mata: __x = st_matrix("__x"); st_matrix("__AD", abs(__x * J(1, `n', 1) - J(`n', 1, 1) * __x'))
		restore
		local opt "edgecov(__AD)"
	}
	local K = rowsof(RB)
	matrix SB = J($NS, `K', .)
	matrix SM = J($NS, `K', .)
	local ok 1
	forvalues s = 1/$NS {
		quietly nwlsm net_`case', dim(2) `opt' clusters(`g') seed(`=2000 + `s'')
		matrix b = e(b)
		matrix T = e(summary)
		forvalues k = 1/`K' {
			matrix SB[`s', `k'] = b[1, `k']
			matrix SM[`s', `k'] = T[`k', 3]
		}
		matrix Q = e(Q)
		matrix CL = e(class)
		mata: st_local("ari", strofreal(lsmt_ari(st_matrix("CL"), st_matrix("TR")), "%7.4f"))
		mata: st_local("ari2", strofreal(lsmt_ari(st_matrix("CL"), st_matrix("LC")), "%7.4f"))
		mata: st_local("rs", strofreal(max(abs(rowsum(st_matrix("Q")) :- 1)), "%12.3e"))
		di as text "case `case' seed `s': ARI with the planted clusters " `ari' "  with latentnet " `ari2' "  max |row sum - 1| " `rs'
		if `ari' < 0.9 | `ari2' < 0.9 | `rs' > 1e-9 local ok 0
		* the classification is the cluster of the largest probability
		mata: st_local("bad", strofreal(lsmt_cls_bad(st_matrix("Q"), st_matrix("CL"))))
		if `bad' > 0 local ok 0
	}
	forvalues k = 1/`K' {
		mata: st_local("avg", strofreal(mean(st_matrix("SB")[., `k'])))
		mata: st_local("se", strofreal(sqrt(mean(st_matrix("SM")[., `k'] :^ 2) / $NS + (st_matrix("RB")[`k', 1] - st_matrix("RB2")[`k', 1]) ^ 2 / 2)))
		local z = abs(`avg' - RB[`k', 1]) / max(`se', 1e-6)
		di as text "case `case' coef `k': nwlsm " %8.4f `avg' "  latentnet " %8.4f RB[`k', 1] "  (posterior sd " %6.4f RB[`k', 2] ")  diff/MCse = " %5.2f `z'
		if `z' > 3.5 local ok 0
	}
	* chains and threads
	quietly nwlsm net_`case', dim(2) `opt' clusters(`g') seed(77) chains(3) threads(2)
	matrix CL = e(class)
	mata: st_local("ari", strofreal(lsmt_ari(st_matrix("CL"), st_matrix("TR")), "%7.4f"))
	di as text "case `case' chains(3): ARI with the planted clusters " `ari'
	if `ari' < 0.9 local ok 0
	if !`ok' {
		di as error "FAIL case `case'"
		exit 9
	}
	di as result "PASS case `case'"
end
lsmt_clu clu45, n(45) dir(0) g(3)
lsmt_clu clu36, n(36) dir(1) g(2) x
display "PASS test_nwlsm_clusters: planted clusters recovered (ARI >= 0.9), coefficients agree with latentnet"
