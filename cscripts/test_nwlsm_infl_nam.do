* Influence model against the linear network autocorrelation (exposure) model with a fixed row-normalised weight matrix: core nwspmatrix + Stata's spregress
* (y on the lagged outcome and the spatially lagged lagged outcome, ivarlag()), which equals least squares with the constructed regressor W ylag.
* Design 1: the true weights are the observed network (y = b0 + b1 ylag + b2 W_A ylag + e, ties generated from latent positions in clusters), b2 = 0.7:
*   both estimates of b2 must have the sign of the truth and agree in magnitude (the latent-distance weights are a smooth version of the observed network).
* Design 2: the true weights are the latent-distance weights of the influence model (5 nearest actors, inverse distance), b2 = 1: the influence model must be
*   about unbiased, the fixed-network model is biased and has a larger root mean square error.
* REPS replications per design (default 20; env NWLSM_REPS).
capture which nwspmatrix
if _rc {
	di "SKIP: test_nwlsm_infl_nam needs the core command nwspmatrix (nwcommands 3.28.0 or later)"
	exit
}
do cscripts/nwlsm_testlib.do
do cscripts/nwlsm_infl_sim.do
local reps : env NWLSM_REPS
if "`reps'" == "" local reps 20
mata:
mata set matastrict off
// design 1 data: a clustered directed network, outcome from the row-normalised observed network
void lsmi_sim_nam(real scalar n, real scalar b2, real scalar seed)
{
	real matrix Z, D, A, W, C
	real colvector ylag, y, cl
	real scalar k
	rseed(seed)
	C = rnormal(4, 2, 0, 3)
	cl = ceil(runiform(n, 1) * 4)
	Z = C[cl, .] + rnormal(n, 2, 0, 1)
	D = lsmt_pd(Z)
	A = (runiform(n, n) :< invlogit(0.5 :- D))
	_diag(A, 0)
	W = A :/ (rowsum(A) :+ (rowsum(A) :== 0))
	ylag = rnormal(n, 1, 0, 1) + 0.8 * rnormal(4, 1, 0, 1)[cl]
	y = 0.5 :+ 0.6 * ylag + b2 * (W * ylag) + rnormal(n, 1, 0, 0.4)
	st_store(., (1..n), A)
	st_store(., n + 1, y)
	st_store(., n + 2, ylag)
	st_store(., n + 3, rnormal(n, 1, 0, 1))
}
end
capture program drop lsmi_fitnam
* least squares / spregress for y on ylag and W ylag; returns the estimate of the W ylag coefficient in r(b)
program lsmi_fitnam, rclass
	syntax, NET(string)
	capture drop _nwid
	nwspmatrix `net', name(Wn) replace silent
	capture quietly spregress y ylag, ml ivarlag(Wn: ylag)
	if _rc {
		di as error "spregress with ivarlag() failed (rc `=_rc')"
		exit 9
	}
	return scalar b = _b[Wn:ylag]
	return scalar rc = 0
end
capture program drop lsmi_nam
program lsmi_nam
	syntax, DESIGN(integer) N(integer) B2(real) REPS(integer) SEED(integer)
	tempname R
	matrix `R' = J(`reps', 3, .)
	forvalues r = 1/`reps' {
		if `design' == 1 {
			clear
			set obs `n'
			forvalues j = 1/`n' {
				quietly gen double __a`j' = .
			}
			quietly gen double y = .
			quietly gen double ylag = .
			quietly gen double x = .
			mata: lsmi_sim_nam(`n', `b2', `=`seed' + `r'')
		}
		else lsmi_sim, n(`n') dir(1) d(2) a0(1.5) b0(0.5) b1(0.6) b2(`b2') b3(0) s2(0.16) g(5) wt(0) seed(`=`seed' + `r'') zsd(2.5)
		nwset __a*, name(cmp) directed replace
		lsmi_fitnam, net(cmp)
		matrix `R'[`r', 1] = r(b)
		quietly nwlsm cmp, dim(2) influence(y ylag) neighbors(5) weight(invdist) seed(`=`seed' + 500 + `r'') nsample(1500) burnin(2500) interval(5) threads(1)
		matrix `R'[`r', 2] = e(b)[1, colsof(e(b))]
		matrix `R'[`r', 3] = (`b2' >= e(summary)[rowsof(e(summary)), 4] & `b2' <= e(summary)[rowsof(e(summary)), 5])
	}
	mata: st_matrix("__m", mean(st_matrix("`R'")))
	mata: st_matrix("__r", (sqrt(mean((st_matrix("`R'")[., 1] :- `b2') :^ 2)), sqrt(mean((st_matrix("`R'")[., 2] :- `b2') :^ 2))))
	mata: st_matrix("__c", correlation(st_matrix("`R'")[., 1..2])[1, 2])
	di as text "DESIGN `design': n=`n', true b2 = `b2', `reps' replications"
	di as text "  mean estimate: fixed network (spregress) " %7.4f __m[1,1] "   influence model " %7.4f __m[1,2] "   (95% interval covers the truth in " %4.2f __m[1,3] " of the replications)"
	di as text "  root mean square error: fixed network " %7.4f __r[1,1] "   influence model " %7.4f __r[1,2] "   correlation of the two estimates " %6.3f __c[1,1]
	if `design' == 1 {
		assert __m[1,1] > 0 & __m[1,2] > 0
		assert abs(__m[1,2] - __m[1,1]) < 0.35
		di as result "PASS design 1: both models find the influence effect with the sign and a similar magnitude"
	}
	else {
		assert abs(__m[1,2] - `b2') < 0.25
		assert __r[1,2] < __r[1,1]
		assert __m[1,3] >= 0.8
		di as result "PASS design 2: the influence model recovers the latent-distance weights, the fixed network model does not"
	}
end
lsmi_nam, design(1) n(100) b2(0.7) reps(`reps') seed(8100)
lsmi_nam, design(2) n(100) b2(1) reps(`reps') seed(8200)
display "PASS test_nwlsm_infl_nam"
