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

## Valued families: details worth knowing
- Plugin blob: for a valued family the y block is values, then an observed mask, then trials (binomial). The plugin protocol was raised; stale plugins are refused by `_nwlsm_libcheck`. Only the macOS plugin is rebuilt locally; the Linux and Windows plugins are stale until CI rebuilds them after a push touching `native/`.
- Normal family: latentnet requires `prior.var` and `prior.var.df`; our defaults 1 and 3 are our own. The MKL objective is divided by the posterior mean dispersion; latentnet's MKL dispersion is probably invalid, so our objective can differ slightly (better or equal on the same expected values).
- latentnet 2.12.0 `simulate()` for UNDIRECTED networks draws every ordered pair and keeps a tie if either draw is a tie (about 1.5 times the edges), so its undirected `gof()` is biased; `estat gof` draws once per pair. The undirected reference in `dev/ref_latentnet.R` (`gof_fixture_u`) is built from latentnet's draws with one draw per dyad. latentnet `gof(model)` fails for latent formulas.
- Mata: `rpoisson(r,c,m)` with a matrix `m` returns rows(m)*r rows, use `rpoisson(1,1,m)`; `import`/`preserve` clear `r()`, stash `r()` matrices first; a ternary with a 1x1 matrix condition gave r(3204), use `if`; a `'` just before `}` on one line breaks do-file parsing; Stata matrix row names cannot contain spaces.
- `estat gof` is deterministic given `seed()` and does not depend on `threads()` (per-draw generator seed + 7919 (k+1)).

## Influence model (0.6.0)
- Identification: with all other actors as neighbours the weighted means m_k of the lagged outcomes are nearly constant (they average over everybody), so the influence effect is collinear with the intercept and its posterior is wide; a moderate `neighbors()` identifies it better. A `neighbors()` much larger than the true number of influential actors biases the effect (paper, section 6.3).
- The tie-model intercept and the position variance are estimated (latentnet prior), the paper fixes the intercept at 2. Distances, and so the weights, are in the units of the tie model with distance coefficient 1.
- `e(b)` has two equations (`network`, `outcome`, both with a `_cons`); the code that reads the tie coefficients uses `e(covariates)` (network names only), not `colnames e(b)`. `e(tuning)` has the network coefficients only.
- `e(Z)` is the minimum Kullback-Leibler fit of the tie model only; the outcome informs the draws (and so `e(D)`, `e(Zmean)`, `e(W)`) but not the MKL optimisation.
- The sets of nearest actors are defined by `d_kj <= G-th smallest distance` of row k; the shortcut updates in `infl_dll()` are valid only for the cases listed there (the mover is in the set and not the boundary and stays inside; or is outside and stays outside); anything else recomputes row k. The sums are rebuilt at the start of every sweep.
- `state_alloc()` is called on a zeroed `State` everywhere (`memset`), including in `stata_call()`; the influence arrays are only allocated with the influence model and freed by `state_free()`.
- Dialog `.dlg` `if` conditions must not use parentheses; the macOS dialog assembly harness is used to check the generated commands.
- `LsmCheckInfl()` (Mata) and plugin codes 5 and 6 exist only for the tests.

## BIC (0.7.0)
- `estat bic` needs the Mata externals of the last fit (like `estat gof`) and `e(mkl_effects)` (only stored when the actor matrices fit into a Stata matrix). `__nwlsm_MZ` holds the MKL positions after centring and rotation; the BIC of the positions without clusters depends on the centring, so keep `lsm_ztrans`.
- The mixture EM is seeded internally (`rseed`, state restored) so it does not change the user's random numbers.
- The worker team must not be used by two callers at once; the chains use one team each, the MKL fit one team (single start only).

- A different `estat bic` with sender and receiver effects than latentnet (442.8 against 399.2 on the monks) is not a BIC-code difference: fed with latentnet's own MKL effects the code reproduces latentnet to 1e-3 (test_nwlsm_bic, bicsamp_sr). The MKL objective has two local optima for weakly identified effects (sd 0.27 or 0.65); latentnet lands in the bad one in 3 of 12 seeds too (`dev/ref/bicsamp_sr_seeds.csv`). nwlsm uses 4 MKL starts by default with effects (deliberate deviation, latentnet uses one), explicit `mklstarts(1)` reproduces the one-start behaviour; the effects part of the BIC is `n log(mean(effect^2))`, so it jumps with the optimum.
