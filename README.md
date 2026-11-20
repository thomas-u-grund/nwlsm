# nwlsm

<!--version-->**Version 1.0.0** (released 20 November 2026)<!--/version-->

Latent space models for networks in Stata, fitted by MCMC with a native C plugin. Part of the nwcommands family; requires the core package
`nwcommands` (3.26.0 or later). Free software, GPL-3. Follows the latent space model of Hoff, Raftery and Handcock (2002) as implemented in
the R package latentnet; a later stage adds the latent space model for social influence of Sweet and Adhikari (2020).

Status: development, not released. Implemented: the Euclidean and bilinear latent space models (`nwlsm netname, dim() edgecov() absdiff() sender receiver sociality bilinear`),
binary and valued (Poisson, binomial, normal: `family() trials()`) undirected and directed networks, random sender, receiver and sociality effects, latent clusters (`clusters()`), parallel chains with Gelman-Rubin diagnostics and threads (`chains()`, `threads()`),
minimum Kullback-Leibler positions, `estat positions | mcmcdiag | clusters | predict | gof | ppcheck | bic`, `nwlsmselect` (number of clusters by BIC as latentnet), `nwlsmplot` and dialogs, validated against latentnet (see `cscripts/`).
The latent space model for social influence of Sweet and Adhikari (2020) is implemented as an option of the same command: `nwlsm netname, influence(outcome lagged_outcome [covariates]) neighbors(#) weight(invdist|expdist)` fits the tie model and a linear model for an actor-level outcome whose
regressors include the distance-weighted lagged outcomes of the other actors, jointly by MCMC; `estat influence` summarises the influence effect, the predictions and who exerts influence. Validated by simulation from the exact model and against the linear network autocorrelation model (see `cscripts/`).
See `dev/DESIGN.md` for the plan.

Reference: Sweet, T. M. and Adhikari, S. 2020. A latent space network model for social influence. Psychometrika 85: 251-274.
