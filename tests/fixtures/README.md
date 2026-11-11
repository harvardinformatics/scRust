# tests/fixtures

Small fixtures committed to the repository, which `cargo test` checks on every change (docs/PLAN.md A1, A4 and A5b): tier-T data exported from the workflow's test data, and R-generated vectors for `rcompat` (stage S0).

- **Size:** about 20 MB in total at most. Larger fixtures go to `data/` (gitignored) and are listed in `validation/manifest.toml`.
- **Layout:** one directory per crate or tool, for example `rcompat/` (R vectors, M1) and `doubletfinder/<object>/` (inputs and stage references for the five test objects, M2).
- **Provenance:** every file is listed in `validation/manifest.toml` with its sha256, size, the script that made it, and the env and profile it was made in. Regenerate fixtures from those recipes; never edit them by hand.
- **Profile:** fixtures that depend on BLAS are made in the portable profile (`OPENBLAS_CORETYPE=Haswell`, one thread), so tests can check them exactly on any x86-64 CPU with AVX2 (A5b).
- **From a crate's tests:** `concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/fixtures")`.

## rcompat/

R 4.5.3's results for the `rcompat` crate (stage S0), made by [`validation/R/rcompat_vectors.R`](../../validation/R/rcompat_vectors.R) and checked by `crates/rcompat/tests/r_vectors.rs`. One tab-separated file per family (`unif`, `random_seed`, `random_seed_after`, `sample`, `rnorm`, `qnorm`, `round`, `pow`, `seq_by`, `seq_len`, `order`, `summary`, `format`); lines starting with `#` name the columns. Doubles are written as their IEEE-754 bit patterns (16 hex digits), so the comparison is bitwise; lists are comma-separated, `-` for an empty one.
