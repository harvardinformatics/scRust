//! The `scrust` binary: `--version`, `fingerprint` (environment fingerprint, M2), and
//! `doubletfinder`, `validate` and `knn` (M3).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn scrust(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_scrust"))
        .args(args)
        .output()
        .expect("failed to run the scrust binary")
}

#[test]
fn version() {
    let out = scrust(&["--version"]);
    assert!(out.status.success());
    assert_eq!(
        String::from_utf8_lossy(&out.stdout),
        format!("scrust {}\n", env!("CARGO_PKG_VERSION"))
    );
}

#[test]
fn doubletfinder_needs_stages_until_sctransform_and_pca_are_ported() {
    let f = fixture("soupx_tenx_test");
    let out = scrust(&[
        "doubletfinder", "--counts", &f.join("input").display().to_string(),
        "--clusters", &f.join("input/seurat_clusters.tsv").display().to_string(),
        "--out", &tmp("no_stages").display().to_string(),
    ]);
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("--from-stages is required"));
}

fn fixture(object: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/doubletfinder").join(object)
}

fn tmp(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli").join(name);
    let _ = std::fs::remove_dir_all(&d);
    d
}

/// Teacher-forced run on a tier-T fixture: the same files as R's reference run, and
/// `scrust validate` finds nothing that differs.
#[test]
fn doubletfinder_reproduces_a_reference_run() {
    let f = fixture("soupx_emptydrops_test");
    let dir = tmp("soupx_emptydrops_test");
    let out = scrust(&[
        "doubletfinder", "--counts", &f.join("input").display().to_string(),
        "--clusters", &f.join("input/seurat_clusters.tsv").display().to_string(),
        "--from-stages", &f.join("stages").display().to_string(),
        "--out", &dir.display().to_string(), "--dump-stages", "light", "--threads", "3",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    for file in ["cells.tsv", "sweep_stats.tsv", "bcmvn.tsv"] {
        assert!(std::fs::read(f.join("reference").join(file)).unwrap() == std::fs::read(dir.join(file)).unwrap(), "{file}");
    }
    let summary = std::fs::read_to_string(dir.join("summary.json")).unwrap();
    assert!(summary.contains("\"pANN_column\": \"pANN_0.25_0.28_30\""));
    let report = dir.join("validate.json");
    let v = scrust(&[
        "validate", "doubletfinder", "--reference", &f.join("reference").display().to_string(),
        "--reference-stages", &f.join("stages").display().to_string(),
        "--candidate", &dir.display().to_string(), "--json", &report.display().to_string(),
    ]);
    assert!(v.status.success(), "{}", String::from_utf8_lossy(&v.stderr));
    assert!(std::fs::read_to_string(&report).unwrap().contains("\"pass\": true"));
    // a second run into the same directory is refused
    let again = scrust(&[
        "doubletfinder", "--counts", &f.join("input").display().to_string(),
        "--clusters", &f.join("input/seurat_clusters.tsv").display().to_string(),
        "--from-stages", &f.join("stages").display().to_string(), "--out", &dir.display().to_string(),
    ]);
    assert_eq!(again.status.code(), Some(2));
}

/// `scrust knn --output count` on the final run's embeddings gives pANN's numerators.
#[test]
fn knn_counts_give_the_final_pann() {
    let f = fixture("decontx_raw_tenx_test");
    let dir = tmp("knn");
    std::fs::create_dir_all(&dir).unwrap();
    let counts = dir.join("counts.npy");
    // 300 real cells, k = round(400 * 0.09) = 36
    let out = scrust(&[
        "knn", "--from-stages", &f.join("stages").display().to_string(), "--context", "final",
        "--queries", "300", "--neighbours", "37", "--output", "count", "--out", &counts.display().to_string(),
        "--threads", "2",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let c = scio::Npy::read(&counts).unwrap();
    let pann = scio::Npy::read(&f.join("stages/S5_pann_final.npy")).unwrap();
    for i in 0..300 {
        assert_eq!((c.get(i, 0) / 36.0).to_bits(), pann.get(i, 0).to_bits(), "cell {i}");
    }
}

#[test]
fn fingerprint_reports_the_pinned_openblas_kernel() {
    // OpenBLAS reads OPENBLAS_CORETYPE when it loads, so set it for the child process: the
    // portable profile's kernel, which every x86-64 CPU with AVX2 can run (docs/PLAN.md A5b).
    let out = Command::new(env!("CARGO_BIN_EXE_scrust"))
        .arg("fingerprint")
        .env("OPENBLAS_CORETYPE", "Haswell")
        .env("OPENBLAS_NUM_THREADS", "1")
        .output()
        .expect("failed to run the scrust binary");
    assert!(out.status.success());
    let json = String::from_utf8_lossy(&out.stdout);
    for line in [
        "\"openblas_corename\": \"Haswell\"",
        "\"openblas_num_threads\": 1",
        "\"OPENBLAS_CORETYPE\": \"Haswell\"",
    ] {
        assert!(json.contains(line), "missing {line} in:\n{json}");
    }
    assert!(json.contains("\"glibc\": \"2."), "no glibc version in:\n{json}");
    assert!(json.contains("libopenblas"), "OpenBLAS library not found in:\n{json}");
}

#[test]
fn unknown_command_is_a_usage_error() {
    let out = scrust(&["no-such-tool"]);
    assert_eq!(out.status.code(), Some(2));
}

/// SHA-256 (FIPS 180-4), for comparing a large output with the digest R recorded.
fn sha256_hex(data: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];
    let mut h: [u32; 8] = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];
    let mut msg = data.to_vec();
    let bits = (data.len() as u64).wrapping_mul(8);
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bits.to_be_bytes());
    for chunk in msg.chunks(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([chunk[4 * i], chunk[4 * i + 1], chunk[4 * i + 2], chunk[4 * i + 3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
        }
        let mut v = h;
        for i in 0..64 {
            let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
            let ch = (v[4] & v[5]) ^ (!v[4] & v[6]);
            let t1 = v[7].wrapping_add(s1).wrapping_add(ch).wrapping_add(K[i]).wrapping_add(w[i]);
            let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
            let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
            let t2 = s0.wrapping_add(maj);
            v = [t1.wrapping_add(t2), v[0], v[1], v[2], v[3].wrapping_add(t1), v[4], v[5], v[6]];
        }
        for (a, b) in h.iter_mut().zip(v) {
            *a = a.wrapping_add(b);
        }
    }
    h.iter().map(|x| format!("{x:08x}")).collect()
}

/// `scrust sctransform` on a tier-T object's counts gives R's SCTransform() outputs byte for
/// byte (validation/R/sctransform_cli_reference.R; scale.data by its SHA-256). Its products are
/// OpenBLAS calls, so it needs the reference's kernel (env.sh's OPENBLAS_CORETYPE=Haswell).
#[test]
fn sctransform_matches_r_byte_for_byte() {
    assert_eq!(sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    if std::env::var("OPENBLAS_CORETYPE").ok().as_deref() != Some("Haswell") {
        eprintln!("skipped: needs OPENBLAS_CORETYPE=Haswell (source env.sh)");
        return;
    }
    let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/sctransform/cli_soupx_tenx_test");
    let out = tmp("sctransform");
    let run = scrust(&[
        "sctransform", "--counts", &fixture("soupx_tenx_test").join("input").display().to_string(),
        "--out", &out.display().to_string(), "--threads", "3",
    ]);
    assert!(run.status.success(), "{}", String::from_utf8_lossy(&run.stderr));
    for f in ["cell_attr.tsv", "gene_attr.tsv", "model_pars.tsv", "model_pars_fit.tsv", "variable_features.tsv"] {
        let got = std::fs::read(out.join(f)).unwrap();
        let want = std::fs::read(reference.join(f)).unwrap();
        assert!(got == want, "{f} differs from R's");
    }
    let digest = std::fs::read_to_string(reference.join("scale_data.npy.sha256")).unwrap();
    assert_eq!(sha256_hex(&std::fs::read(out.join("scale_data.npy")).unwrap()), digest.trim(), "scale_data.npy differs from R's");
    assert!(out.join("summary.json").is_file());
}
