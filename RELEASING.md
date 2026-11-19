# Releasing nwlsm

The master release procedure for all packages of the family is `RELEASING.md` of the core repository (nwcommands). This file keeps only what is specific to nwlsm.

- The public repository is https://github.com/thomas-u-grund/nwlsm (branch `main`). The CI commits the plugins only after a push to `main` that touches `native/` or the workflow file; a manual run builds but does not commit. Tagging and pushing are the owner's decision.
- Versions: `VERSION` + `python3 tools/set_version.py` / `python3 tools/set_version.py --check`; semantic versioning.
- Build the Mata library with Stata 17: `StataBE -e do lib/build.do` (Stata 17 only; the build stops under a newer Stata unless `NW_ALLOW_NEWER_STATA=1`).
- Build the plugin: `cd native && make macos` (macOS), `make unix` (Linux); Windows by `cl /LD /O2 /DSYSTEM=STWIN32 /I spi lsm_sim.c spi\stplugin.c /Fe:..\lib\plugins\windows\lsm_sim_windows.plugin`. After the first push the workflow `.github/workflows/build-plugins.yml` builds all three and commits the binaries.
- Tests: `NWCORE=<core checkout> cscripts/run_tests.sh`. The slow tests are `test_nwlsm_exact.do` (about 80 s) and `test_nwlsm_karate.do` (about 2 minutes). The fixtures of `test_nwlsm_changestats.do` are made by `dev/gen_ref_definitions.py`.
- Integration into the core (`nwinstall`, `nwversion`, website) is a core release and is not part of this repository.
