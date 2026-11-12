//! Records the conda env's `lib/` as the rpath of this package's binaries and tests, so they find
//! the OpenBLAS that rcompat links (`-lblas`, see `crates/rcompat/build.rs`) at run time.

fn main() {
    println!("cargo:rerun-if-env-changed=CONDA_PREFIX");
    if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
        if !prefix.is_empty() {
            println!("cargo:rustc-link-arg=-Wl,-rpath,{prefix}/lib");
        }
    }
}
