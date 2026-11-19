# nwlsm

<!--version-->**Version 0.3.0** (released 19 November 2026)<!--/version-->

Latent space models for networks in Stata, fitted by MCMC with a native C plugin. Part of the nwcommands family; requires the core package
`nwcommands` (3.26.0 or later). Free software, GPL-3. Follows the latent space model of Hoff, Raftery and Handcock (2002) as implemented in
the R package latentnet; a later stage adds the latent space model for social influence of Sweet and Adhikari (2020).

Status: development, not released. Implemented: the Euclidean and bilinear latent space models (`nwlsm netname, dim() edgecov() absdiff() sender receiver sociality bilinear`),
binary undirected and directed networks, random sender, receiver and sociality effects, minimum Kullback-Leibler positions, `estat positions | mcmcdiag | predict`
and `nwlsmplot`, validated against latentnet (see `cscripts/`). Not yet implemented: latent clusters, valued networks. See `dev/DESIGN.md` for the plan.
