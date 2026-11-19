# TODO nwlsm
See dev/DESIGN.md.
- Influence model, not yet done: the hierarchical version for several networks with random coefficients (Sweet and Adhikari 2020, section 5), the paper's susceptibility indicator (actors far from every other actor are not influenced), several outcome waves, a fixed intercept of the tie model (the paper fixes it at 2), binary or count outcomes (a contagion version, not a published method).
- Persistent worker pool for the per-dyad passes (threads are created for every pass now; a pool would make threads pay from smaller networks on).
- Run the tests on Linux and Windows; read the CI plugins (rebuild after pushing `native/`).
- Help file examples with a real data set; website and core integration are a core release.
- Cluster number selection (BIC-type criterion as in latentnet): needs the marginal likelihood of the mixture models (latentnet's BICmcmc and its approximation by the posterior means); not cheap here, so only the log likelihoods of fits with different clusters() can be compared by the user.
- `estat gof` for valued networks uses the value distribution only; geodesic statistics of a thresholded network could be added.
