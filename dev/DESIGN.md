# nwlsm design (internal; not shipped)

Owner: Thomas U. Grund (sole author). GPL-3. Part of the nwcommands family; depends only on the core. Native C plugin only (no Mata fallback).

## Status (2026-11-19)
Stage 1, the sender/receiver/sociality and bilinear parts and the latent clusters of stage 2 are implemented (0.4.0, with chains and threads); valued networks are open.

## Scope in stages
1. **Base model (latentnet "euclidean")**: binary one-mode networks (undirected and directed), logit P(y_ij = 1) = beta0 + sum_k beta_k x_ijk - |z_i - z_j|, z_i in R^d. Posterior by MCMC (random-walk Metropolis-Hastings on beta and z; Gibbs for the latent-position variance). Output: posterior means and intervals of beta, posterior sample, positions after Procrustes alignment, and the minimum-Kullback-Leibler (MKL) positions, as latentnet. Edge covariates from nwcommands networks.
2. Extensions that latentnet offers, in order of value: sender/receiver random effects (`sender`, `receiver`, `rsociality`), bilinear distance (`bilinear`), latent clusters (`G=`), valued networks (Poisson, ordinal). Each only after the previous one reproduces latentnet.
3. **Latent space model for social influence** (Sweet and Adhikari 2020, Psychometrika 85(2): 251-274): a Gaussian outcome at time t regressed on its lag and on the distance-weighted mean of the lagged outcomes of the G nearest actors in the latent space, estimated jointly with the latent space model of the ties. No reference software known (the authors wrote their own sampler); validate by simulation and by the numbers in the paper. Read the paper in full before implementing (we have read the model section and the first simulation example only).
4. Possible later: a cross-sectional binary-outcome version for the contagion study (not a published method; only if the paper needs it).

## Reference for stage 1
R package latentnet 2.12.0 (installed in the owner's R library; GPL). Defaults of `ergmm` (checked on 2026-11-19 with data(sampson), d = 2, no clusters, `fit$prior`): beta ~ N(0, var = 9) with `adjust.beta.var = TRUE` (the prior variance of a covariate coefficient is divided by the mean square of that covariate); Z variance prior (unclustered model): z_i ~ N(0, sigma2 I_d), sigma2 scale-inverse-chi-squared with scale Z.var = 0.125 * n^(2/d) and degrees of freedom Z.var.df = sqrt(n) (checked for n = 18, 20, 40, 80 and d = 2, 3 in dev/ref/prior_formulas.csv); Z.mean.var = 2 * Z.var (clusters only); sampling defaults: burn-in 10,000, 4,000 draws (sample.size), interval 10, target acceptance 0.234 with pilot runs tuning Z.delta = 0.6.
Compare: posterior means of beta, the posterior mean of the log likelihood, the matrix of posterior mean distances (invariant to rotation, reflection and translation, so it can be compared directly; the positions cannot), and the MKL positions after Procrustes. Because both are MCMC, agreement is statistical: use several seeds and compare with Monte Carlo error.

## Validation plan
- Likelihood: log-likelihood for fixed (beta, Z) equals an independent calculation in Stata/Mata and in R.
- Single-block Metropolis updates: acceptance ratios equal the full recomputation of the log posterior difference (test with random states).
- Recovery: simulate networks from known (beta, Z), check coverage of beta and recovery of the distance matrix.
- Agreement with latentnet on: Sampson (data(sampson), directed), a random undirected network with a node covariate, a larger network (n = 100).
- Speed: report time against latentnet on n = 50, 100, 200.

## Conventions to follow (from the core RELEASING.md / nwalaam)
Version tooling (tools/set_version.py, check_version.py), VERSION file with the date, CHANGELOG entry per version, native/Makefile plus the workflow in .github/workflows (Linux and Windows plugins are built by CI on a push touching native/), plugin protocol number checked by the .ado (never shown to users), tests in cscripts/ run with the scratch-PLUS pattern of CLAUDE.md, user-facing text never cites internal files. No push, tag or release without the owner's go.

## Timing (2026-11-19, macOS, Stata 17 BE vs R 4.x latentnet 2.12.0, default settings: 10,000 burn-in, 4,000 draws, interval 10; directed network, d = 2, sender and receiver effects, MKL included)
| n | latentnet | nwlsm |
|---|---|---|
| 50 | 7.1 s | 4.4 s |
| 100 | 31.4 s | 19.0 s |
| 200 | 117.8 s | 76.8 s |
The run is bound by exp/log1p per dyad (about 7.7 ns per evaluation, 1e10 evaluations at n = 200); -O3 -fno-math-errno did not help. The n updates of the sender effects (and then of the receiver effects) are conditionally independent and could run in threads.
