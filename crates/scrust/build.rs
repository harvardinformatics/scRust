//! Links the `scrust` binary against the conda env's OpenBLAS (docs/PLAN.md A5 and A5b), so
//! `scrust fingerprint` reports the kernel the same OpenBLAS build chooses in R.
//!
//! The conda packages ship no unversioned `libopenblas.so`, so link with `-lblas`
//! (`libblas.so -> libopenblasp-r0.3.34.so`, the library R uses), not `-lopenblas`
//! (envs/rust-dev.yml). The env's `lib/` is taken from `CONDA_PREFIX` (set by `source env.sh`)
//! and recorded as the binary's rpath, so the binary finds that same library at run time.

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
}
