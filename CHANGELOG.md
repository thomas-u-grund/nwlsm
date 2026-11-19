# Changelog

## [0.2.0] - 2026-11-19
- First estimation command: `nwlsm netname, dim() edgecov() absdiff() nsample() burnin() interval() seed()` fits the Euclidean latent space model of latentnet (binary one-mode networks, undirected and directed) by MCMC in the native plugin `lsm_sim` (random-walk Metropolis-Hastings for the positions, one actor at a time, and for the coefficients; Gibbs update of the position variance; latentnet's default priors including the covariate-adjusted coefficient variance; step sizes tuned in pilot runs toward an acceptance rate of 0.234). Replay, `e()` results, posterior mean distance matrix `e(D)`, Procrustes-aligned mean positions `e(Z)`, `savedraws()`.
- Minimum Kullback-Leibler (MKL) positions of latentnet are NOT implemented yet; positions are the mean of the Procrustes-aligned draws.
- Not implemented in this version: sender/receiver effects, bilinear distance, clusters, valued networks, `estat`, plots.
- Run so far (Stata 17 BE on macOS only): `test_nwlsm_check` (log likelihood and one-actor update agree with a direct recomputation in 60 random states), `test_nwlsm_latentnet` (Sampson d = 2, random n = 40 with d = 2 and d = 3 against the latentnet reference), `test_nwlsm_recovery`, `test_nwlsm_errors`. Not run: Stata 18/19, Windows, Linux (the Linux and Windows plugins are built by the workflow in `.github/workflows`).

## [0.1.0] - 2026-11-19
- Repository created: version tooling, build workflow, help stub, design notes. No estimation command yet.
