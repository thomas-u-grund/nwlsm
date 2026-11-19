# GOTCHA - known pitfalls of nwlsm

## The MKL positions are not unique, and latentnet's do not have to agree with ours
The objective of the minimum Kullback-Leibler fit (likelihood of the posterior mean tie probabilities) is flat or multimodal in the positions, most of all with sender/receiver effects and in sparse networks. latentnet starts a BFGS from the single closest draw and ends in a local optimum; for Sampson with sender and receiver effects its positions correlate only 0.23 with its own posterior mean distances (checked 2026-11-19; its objective differs from the best of ten random starts by 0.1 log-likelihood units). So the tests compare the objective (recomputed from the returned parameters, not below the objective at the posterior means) for models with effects, and the pairwise distances against latentnet only for the plain Euclidean fixtures, with a threshold of 0.90. `mklstarts(1)` (the default since 0.4.0) is latentnet's procedure; `mklstarts(#)` starts from the best # draws and keeps the best optimum.

## The MKL fit needs a ridge
Without a penalty the optimiser ran away along a flat direction (one actor at x = 3e6) on a random network. The plugin adds 1/1000 of the posterior precision (1 / mean s2, 1 / mean effect variance) as a ridge on the positions and effects; `e(mkl_obj)` is the unpenalised likelihood at the optimum. latentnet has none.

## latentnet's centring of the MKL fit
`.scale.ergmm.model` moves the means of the effects into the intercept; for the sociality effect it uses the factor `1 + directed`, i.e. 1 for an undirected network although sociality enters both ends of a tie. nwlsm uses 2 (the intercept shift that leaves the linear predictor unchanged).

## Bilinear model: do not centre
The inner product is not translation-invariant. `lsm_ztrans()` and the Procrustes alignment centre the positions only for the Euclidean model; the bilinear positions are only rotated or reflected.

## Plugin output layout
The result rows of the chain call are laid out in `native/lsm_sim.c` (header of 2K + 8 rows, then draws, then the MKL block) and read in `LsmFit()` in `unw_lsm.do`. Change both together; a mismatch shows up as missing values ("the sampler returned missing values").

## Stata
- `nw_labs_string()` (the `labs` local of `_nwsyntax`) separates the actor names by commas; `matrix rownames` needs spaces.
- A colour style such as `mcolor(p1%70)` drew black markers here (batch, `-e`); use named colours with opacity.
- `st_addvar("double", "_cons")` is an invalid name: `LsmSaveDraws()` renames it `beta_cons`.
- `estat mcmcdiag` reads the draws from the Mata externals of the last fit; they are gone after `clear all` or `mata: mata clear`.
- The `weak ridge` and the multiple starts are deviations from latentnet and are documented in the help file.

## Threads: why results do not depend on threads()
Work items that may run on different threads are independent (chains, blocks of draws for the Kullback-Leibler divergences, MKL starts) or are cut into 16 fixed blocks of rows whose partial sums are added in block order (`beta_pass`, `mkl_fg`); never let a thread add to a shared sum. Chain k has its own generator (seed + 1000003 k) and its own state. No `SF_*` call from a worker thread. `test_nwlsm_threads` compares everything bit for bit.
Threads are created per pass (no pool), so that inner threads pay only from about n = 100 (`n*n/10000` threads at most); the updates of the positions and effects are Gauss-Seidel sweeps of O(n) steps and cannot be split without changing the sampler.

## Chains: positions cannot be pooled directly
Draws from different chains are rotated differently. Everything pooled across chains is rotation-invariant (coefficients, log likelihood, variances, distance and tie probability means, cluster labels after correction); the positions are aligned to the MKL positions (`Zmean`) draw by draw. R-hat is therefore computed on invariant quantities only.

## Cluster labels
The chain's labels switch freely in a cluster model. `lsm_relabel()` permutes the labels of every draw (all G! permutations, G <= 6) to maximise agreement with the reference memberships. `PI` also reorders the variances, the mixing proportions and the cluster means. Mata colon operators do not broadcast a column against a row (`c :== r` gives r(3200)); expand with `J()`.
