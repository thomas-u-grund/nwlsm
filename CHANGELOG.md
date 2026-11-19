# Changelog

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
