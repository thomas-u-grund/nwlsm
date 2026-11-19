# nwlsm design (internal; not shipped)

Owner: Thomas U. Grund (sole author). GPL-3. Part of the nwcommands family; depends only on the core. Native C plugin only (no Mata fallback).

## Scope in stages
1. **Base model (latentnet "euclidean")**: binary one-mode networks (undirected and directed), logit P(y_ij = 1) = beta0 + sum_k beta_k x_ijk - |z_i - z_j|, z_i in R^d. Posterior by MCMC (random-walk Metropolis-Hastings on beta and z; Gibbs for the latent-position variance). Output: posterior means and intervals of beta, posterior sample, positions after Procrustes alignment, and the minimum-Kullback-Leibler (MKL) positions, as latentnet. Edge covariates from nwcommands networks.
2. Extensions that latentnet offers, in order of value: sender/receiver random effects (`sender`, `receiver`, `rsociality`), bilinear distance (`bilinear`), latent clusters (`G=`), valued networks (Poisson, ordinal). Each only after the previous one reproduces latentnet.
3. **Latent space model for social influence** (Sweet and Adhikari 2020, Psychometrika 85(2): 251-274): a Gaussian outcome at time t regressed on its lag and on the distance-weighted mean of the lagged outcomes of the G nearest actors in the latent space, estimated jointly with the latent space model of the ties. No reference software known (the authors wrote their own sampler); validate by simulation and by the numbers in the paper. Read the paper in full before implementing (we have read the model section and the first simulation example only).
4. Possible later: a cross-sectional binary-outcome version for the contagion study (not a published method; only if the paper needs it).

## Reference for stage 1
R package latentnet 2.12.0 (installed in the owner's R library; GPL). Defaults of `ergmm` (checked on 2026-11-19 with data(sampson), d = 2, no clusters, `fit$prior`): beta ~ N(0, var = 9) with `adjust.beta.var = TRUE` (the prior variance of a covariate coefficient is divided by the mean square of that covariate); Z variance prior: scale-inverse-chi-squared with scale Z.var = 2.25 and df Z.var.df = 4.24 for n = 18 (Z.var.mul = 0.125, Z.var.df.mul = 1: apparently Z.var = 0.125 * n and df = sqrt(n); confirm with other n and d with dev/ref_latentnet.R before relying on it); sampling defaults: burn-in 10,000, 4,000 draws (sample.size), interval 10, target acceptance 0.234 with pilot runs tuning Z.delta = 0.6.
Compare: posterior means of beta, the posterior mean of the log likelihood, the matrix of posterior mean distances (invariant to rotation, reflection and translation, so it can be compared directly; the positions cannot), and the MKL positions after Procrustes. Because both are MCMC, agreement is statistical: use several seeds and compare with Monte Carlo error.

## Validation plan
- Likelihood: log-likelihood for fixed (beta, Z) equals an independent calculation in Stata/Mata and in R.
- Single-block Metropolis updates: acceptance ratios equal the full recomputation of the log posterior difference (test with random states).
- Recovery: simulate networks from known (beta, Z), check coverage of beta and recovery of the distance matrix.
- Agreement with latentnet on: Sampson (data(sampson), directed), a random undirected network with a node covariate, a larger network (n = 100).
- Speed: report time against latentnet on n = 50, 100, 200.

## Conventions to follow (from the core RELEASING.md / nwalaam)
Version tooling (tools/set_version.py, check_version.py), VERSION file with the date, CHANGELOG entry per version, native/Makefile plus the workflow in .github/workflows (Linux and Windows plugins are built by CI on a push touching native/), plugin protocol number checked by the .ado (never shown to users), tests in cscripts/ run with the scratch-PLUS pattern of CLAUDE.md, user-facing text never cites internal files. No push, tag or release without the owner's go.
