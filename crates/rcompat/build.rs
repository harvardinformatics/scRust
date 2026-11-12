//! Links rcompat against the conda env's OpenBLAS (docs/PLAN.md A5 and A5b): R's matrix products
//! (`tcrossprod` in sctransform's residuals, `%*%` and `svd` in irlba), Armadillo's `gemv` in
//! glmGamPoi's Cox-Reid term and the PORT library's `ddot`/`dnrm2` all go through this library,
//! so the ports call it too.
//!
//! The conda packages ship no unversioned `libopenblas.so`, so link with `-lblas`
//! (`libblas.so -> libopenblasp-r0.3.34.so`, the library R uses), as `crates/scrust/build.rs` does.
//! The search path and the library propagate to every crate that depends on rcompat; the rpath is
//! a link argument, which applies to this package's own test binaries only (dependents set their
//! own, see `crates/scrust/build.rs`).

fn main() {
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");
    match std::env::var("CONDA_PREFIX") {
        Ok(prefix) if !prefix.is_empty() => {
            let lib = format!("{prefix}/lib");
            println!("cargo:rustc-link-search=native={lib}");
            println!("cargo:rustc-link-arg=-Wl,-rpath,{lib}");
        }
        _ => println!(
            "cargo:warning=CONDA_PREFIX is not set (source env.sh first); \
             linking -lblas from the system library path"
        ),
    }
    println!("cargo:rustc-link-lib=dylib=blas");
    // LAPACK (`dgesdd`, M5): `liblapack.so` is the same OpenBLAS file, whose bundled reference
    // LAPACK R's `La_library()` reports; linked by that name for what it provides
    println!("cargo:rustc-link-lib=dylib=lapack");
}
