# tests/fixtures

Small fixtures committed to the repository, which `cargo test` checks on every change (docs/PLAN.md A1, A4 and A5b): tier-T data exported from the workflow's test data, and R-generated vectors for `rcompat` (stage S0).

- **Size:** about 20 MB in total at most. Larger fixtures go to `data/` (gitignored) and are listed in `validation/manifest.toml`.
- **Layout:** one directory per crate or tool, for example `rcompat/` (R vectors, M1) and `doubletfinder/<object>/` (inputs and stage references for the five test objects, M2).
- **Provenance:** every file is listed in `validation/manifest.toml` with its sha256, size, the script that made it, and the env and profile it was made in. Regenerate fixtures from those recipes; never edit them by hand.
- **Profile:** fixtures that depend on BLAS are made in the portable profile (`OPENBLAS_CORETYPE=Haswell`, one thread), so tests can check them exactly on any x86-64 CPU with AVX2 (A5b).
- **From a crate's tests:** `concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures")`.
