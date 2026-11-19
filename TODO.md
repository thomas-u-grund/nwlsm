# TODO nwlsm
See dev/DESIGN.md.
- Valued networks (latentnet family: Poisson, binomial with trials, Gaussian), not started. Plan: a family code and `trials` matrix in the blob; `llterm()` per family (Poisson: y eta - exp(eta) - lgamma(y+1); binomial: y eta - trials log(1+exp(eta)) + lchoose; Gaussian: -(y-eta)^2/(2 s2) with an inverse-gamma Gibbs step for the dispersion, latentnet's `dispersion` prior); the MKL objective and `EYsum` need the family's mean function; nwlsm must accept valued networks only with `family()`. Stage 3 (social influence).
- Persistent worker pool for the per-dyad passes (threads are created for every pass now; a pool would make threads pay from smaller networks on).
- Run the tests under Stata 18 and 19 and on Linux and Windows; read the CI plugins.
- Help file examples with a real data set; website and core integration are a core release.
- A goodness-of-fit command; cluster number selection (BIC-type criterion as in latentnet).
