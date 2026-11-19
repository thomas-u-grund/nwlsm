# Changelog

## [0.4.0] - 2026-11-19
- `mklstarts()` now defaults to 1: the minimum Kullback-Leibler fit starts from the single closest draw, as in latentnet; several starts remain an option. The very weak ridge of the MKL fit is kept.
- Threads and chains, built into the plugin. `chains(#)` runs independent chains in parallel and pools their draws (`nsample()` is the total over all chains; each chain has its own burn-in and tuning; chain k has the seed `seed()` + 1000003 k; chains after the first start from over-dispersed positions). `estat mcmcdiag` shows the Gelman-Rubin statistic (also `e(rhat)`), effective sample sizes summed over the chains and the largest Geweke z. `threads(#)` (default: one per core) distributes the chains over threads; spare threads work inside a chain on the passes over all dyads (coefficient updates) and on the minimum Kullback-Leibler fit in networks of about 100 actors or more. Results for a given seed and `chains()` are identical for every `threads()` (`test_nwlsm_threads`). pthreads on macOS and Linux, `_beginthreadex` on Windows.
- Latent clusters: `clusters(G)` (G = 2 to 6), latentnet's mixture of spherical normals for the positions with its priors; memberships, cluster means, variances and mixing proportions are sampled by Gibbs steps; labels are corrected over the draws. Results `e(Q)`, `e(class)`, `e(cmeans)`, `e(lambda)`, one variance per cluster; `estat clusters`, `nwlsmplot, clusters`.
- Dialogs: `nwlsm.dlg`, `nwlsmplot.dlg`, `nwlsm_estat.dlg`.
- `savedraws()` also saves the variance of each cluster and the mixing proportions; `LsmSample()` returns `"class"`, `"lambda"` and `"chains"`.
- Tests: `test_nwlsm_clusters` (planted clusters, ARI against truth and latentnet, coefficients against latentnet), `test_nwlsm_threads`; the earlier tests were rerun with the new defaults. Run so far: Stata 17 BE on macOS only.
- Not implemented: valued networks (Poisson, binomial, Gaussian).

## [0.3.0] - 2026-11-19
- Random effects: `sender`, `receiver` (directed networks) and `sociality` (undirected and directed) with latentnet's priors (normal effects, scaled inverse chi-squared variance with scale 1 and 3 degrees of freedom; `sendervar()`, `senderdf()`, ... change them), random-walk Metropolis-Hastings updates of every effect with tuned step sizes and Gibbs updates of the variances, all in the native plugin. Results: `e(sender)`, `e(receiver)`, `e(sociality)`, `e(vsummary)`.
- Minimum Kullback-Leibler (MKL) positions are now the reported positions `e(Z)` (with `e(mkl_b)`, `e(mkl_effects)`), found in the plugin by a limited-memory BFGS maximisation of the likelihood of the posterior mean tie probabilities from the `mklstarts()` (default 5) draws closest to them. Deviations from latentnet: several starts, and a very weak ridge (1/1000 of the posterior precision) on positions and effects against runaway along flat directions. The mean of the draws aligned to the MKL positions is `e(Zmean)`; the posterior mean tie probabilities are `e(P)`.
- Bilinear latent term (`bilinear`, latentnet's `bilinear`).
- Post-estimation: `estat positions`, `estat mcmcdiag` (effective sample sizes and Geweke z-scores), `estat predict` (posterior mean tie probabilities as a valued network), and the plot command `nwlsmplot`.
- `savedraws()` now also saves the variance draws of the effects, and works with the intercept (it failed on the name `_cons` in 0.2.0).
- Tests run so far (Stata 17 BE on macOS only): `test_nwlsm_check` (likelihood, one-actor and one-effect updates, MKL gradient; 90 random states), `test_nwlsm_latentnet` (Sampson, random n = 40 with d = 2 and 3), `test_nwlsm_effects` (five models with sender, receiver, sociality, node covariate and bilinear terms against two latentnet chains each), `test_nwlsm_post`, `test_nwlsm_recovery`, `test_nwlsm_recovery_effects`, `test_nwlsm_errors`. Not run: Stata 18/19, Windows, Linux.
- Not implemented: latent clusters, valued networks.

## [0.2.0] - 2026-11-19
- First estimation command: `nwlsm netname, dim() edgecov() absdiff() nsample() burnin() interval() seed()` fits the Euclidean latent space model of latentnet (binary one-mode networks, undirected and directed) by MCMC in the native plugin `lsm_sim` (random-walk Metropolis-Hastings for the positions, one actor at a time, and for the coefficients; Gibbs update of the position variance; latentnet's default priors including the covariate-adjusted coefficient variance; step sizes tuned in pilot runs toward an acceptance rate of 0.234). Replay, `e()` results, posterior mean distance matrix `e(D)`, Procrustes-aligned mean positions `e(Z)`, `savedraws()`.
- Minimum Kullback-Leibler (MKL) positions of latentnet are NOT implemented yet; positions are the mean of the Procrustes-aligned draws.
- Not implemented in this version: sender/receiver effects, bilinear distance, clusters, valued networks, `estat`, plots.
- Run so far (Stata 17 BE on macOS only): `test_nwlsm_check` (log likelihood and one-actor update agree with a direct recomputation in 60 random states), `test_nwlsm_latentnet` (Sampson d = 2, random n = 40 with d = 2 and d = 3 against the latentnet reference), `test_nwlsm_recovery`, `test_nwlsm_errors`. Not run: Stata 18/19, Windows, Linux (the Linux and Windows plugins are built by the workflow in `.github/workflows`).

## [0.1.0] - 2026-11-19
- Repository created: version tooling, build workflow, help stub, design notes. No estimation command yet.
