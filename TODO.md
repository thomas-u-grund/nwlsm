# TODO nwlsm
See dev/DESIGN.md.
- Stage 3 (social influence, Sweet and Adhikari 2020).
- Persistent worker pool for the per-dyad passes (threads are created for every pass now; a pool would make threads pay from smaller networks on).
- Run the tests on Linux and Windows; read the CI plugins (rebuild after pushing `native/`).
- Help file examples with a real data set; website and core integration are a core release.
- Cluster number selection (BIC-type criterion as in latentnet).
- `estat gof` for valued networks uses the value distribution only; geodesic statistics of a thresholded network could be added.
