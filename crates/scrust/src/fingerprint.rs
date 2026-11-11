//! `scrust fingerprint`: the environment fingerprint of docs/PLAN.md A5b, mechanism 3.
//!
//! Every reference run and every Rust run records the machine and library facts that decide
//! whether results can be compared bit for bit: the node and SLURM job, the CPU model and the
//! instruction-set flags OpenBLAS's kernels need, the OpenBLAS kernel ("core") and thread count
//! actually in use, the kernel-selection variables, and the glibc version. The validation
//! harness runs exact checks only when the BLAS-relevant fields match the reference's.
//!
//! The OpenBLAS facts come from the library itself (`openblas_get_corename`,
//! `openblas_get_config`, `openblas_get_num_threads`), linked as `-lblas` from the conda env
//! (build.rs). R's env uses a byte-identical copy of the same build (envs/rust-dev.yml), so on
//! a given node and with the same `OPENBLAS_*` variables both choose the same kernel; R
//! reference runs call this subcommand for their own fingerprint.

use std::ffi::CStr;
use std::fmt::Write as _;
use std::os::raw::{c_char, c_int};

extern "C" {
    fn openblas_get_corename() -> *const c_char;
    fn openblas_get_config() -> *const c_char;
    fn openblas_get_num_threads() -> c_int;
    fn gnu_get_libc_version() -> *const c_char;
}

/// Environment variables recorded verbatim (`null` when unset).
const ENV_VARS: &[&str] = &[
    "OPENBLAS_CORETYPE",
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "SLURM_JOB_ID",
    "SLURM_JOB_PARTITION",
    "SLURMD_NODENAME",
    "SLURM_CPUS_PER_TASK",
];

/// CPU flags that decide which OpenBLAS kernels can run (Haswell needs avx2 and fma; the
/// SkylakeX, Cooperlake and SapphireRapids kernels need AVX-512).
const CPU_FLAGS: &[&str] = &["avx2", "fma", "avx512f", "avx512bw", "avx512_bf16", "amx_tile"];

fn c_string(p: *const c_char) -> Option<String> {
    if p.is_null() {
        None
    } else {
        // SAFETY: the OpenBLAS and glibc functions return pointers to static, NUL-terminated
        // strings that live as long as the process.
        Some(unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned())
    }
}

fn read_trimmed(path: &str) -> Option<String> {
    std::fs::read_to_string(path)
        .ok()
        .map(|s| s.trim().to_string())
}

/// The first value of `key` in /proc/cpuinfo.
fn cpuinfo_field(cpuinfo: &str, key: &str) -> Option<String> {
    cpuinfo.lines().find_map(|line| {
        let (k, v) = line.split_once(':')?;
        (k.trim() == key).then(|| v.trim().to_string())
    })
}

/// The file the OpenBLAS library was loaded from, as mapped into this process.
fn openblas_library() -> Option<String> {
    let maps = std::fs::read_to_string("/proc/self/maps").ok()?;
    maps.lines().find_map(|line| {
        let path = line.split_whitespace().nth(5)?;
        let file = path.rsplit('/').next()?;
        (file.starts_with("libopenblas") || file.starts_with("libblas")).then(|| path.to_string())
    })
}

fn json_string(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn json_opt(out: &mut String, s: Option<&str>) {
    match s {
        Some(s) => json_string(out, s),
        None => out.push_str("null"),
    }
}

/// The fingerprint as a JSON object (keys in a fixed order, one per line).
pub fn json() -> String {
    // SAFETY: plain C calls with no arguments; see `c_string` for the returned pointers.
    let (corename, config, threads, glibc) = unsafe {
        (
            c_string(openblas_get_corename()),
            c_string(openblas_get_config()),
            openblas_get_num_threads(),
            c_string(gnu_get_libc_version()),
        )
    };
    let cpuinfo = std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
    let flag_line = cpuinfo_field(&cpuinfo, "flags").unwrap_or_default();
    let flags: Vec<&str> = CPU_FLAGS
        .iter()
        .copied()
        .filter(|want| flag_line.split_whitespace().any(|f| f == *want))
        .collect();

    let mut out = String::from("{\n");
    let mut field = |key: &str, value: &dyn Fn(&mut String)| {
        if !out.ends_with("{\n") {
            out.push_str(",\n");
        }
        out.push_str("  ");
        json_string(&mut out, key);
        out.push_str(": ");
        value(&mut out);
    };
    field("scrust_version", &|o| json_string(o, env!("CARGO_PKG_VERSION")));
    field("host", &|o| {
        json_opt(o, read_trimmed("/proc/sys/kernel/hostname").as_deref())
    });
    field("kernel", &|o| {
        json_opt(o, read_trimmed("/proc/sys/kernel/osrelease").as_deref())
    });
    field("cpu_model", &|o| {
        json_opt(o, cpuinfo_field(&cpuinfo, "model name").as_deref())
    });
    field("cpu_flags", &|o| {
        o.push('[');
        for (i, f) in flags.iter().enumerate() {
            if i > 0 {
                o.push_str(", ");
            }
            json_string(o, f);
        }
        o.push(']');
    });
    field("glibc", &|o| json_opt(o, glibc.as_deref()));
    field("openblas_corename", &|o| json_opt(o, corename.as_deref()));
    field("openblas_config", &|o| json_opt(o, config.as_deref()));
    field("openblas_num_threads", &|o| {
        let _ = write!(o, "{threads}");
    });
    field("openblas_library", &|o| json_opt(o, openblas_library().as_deref()));
    for var in ENV_VARS {
        field(var, &|o| json_opt(o, std::env::var(var).ok().as_deref()));
    }
    out.push_str("\n}\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpuinfo_fields() {
        let text = "processor\t: 0\nmodel name\t: Some CPU @ 2.0GHz\nflags\t\t: fpu avx2 fma\n";
        assert_eq!(cpuinfo_field(text, "model name").as_deref(), Some("Some CPU @ 2.0GHz"));
        assert_eq!(cpuinfo_field(text, "flags").as_deref(), Some("fpu avx2 fma"));
        assert_eq!(cpuinfo_field(text, "missing"), None);
    }

    #[test]
    fn json_escapes() {
        let mut s = String::new();
        json_string(&mut s, "a\"b\\c\nd\u{1}");
        assert_eq!(s, "\"a\\\"b\\\\c\\nd\\u0001\"");
    }
}
