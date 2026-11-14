//! Constructed edge cases (docs/PLAN.md C7), against R (`tests/fixtures/clustering/edge/`, made by
//! `validation/R/cluster_references.R --mode edge` and `validation/slurm/fixtures_cluster.sh`):
//! isolated cells, which tie at 0 with every cluster so that R's draw places them (ties of 2, and
//! of 10, which `sample()` resolves to the 9th); singletons whose earlier merges decide the later
//! ones; networks of 63, 64, 65, 128, 256 and 1,024 cells (`nextInt`'s power-of-two branch); one
//! and two cells, disconnected components; equal gains (a ring, a lattice) and equal-sized
//! clusters; labels that become text; every cluster a singleton (R stops); negative and extreme
//! seeds; several resolutions, in any order; all three algorithms; both modularity functions.
//!
//! For each case and setting the test runs FindClusters, recomputes every file R wrote (the labels
//! and levels, the workflow's files, the K7 stage arrays, the RNG checkpoints, the instrumented
//! optimizer's starts, partitions and ids) from Rust's result in R's formats, and compares their
//! sha256 with R's; it compares what FindClusters printed and its messages line by line, and where
//! R stopped, the error.

use cluster::{find_clusters, ClusterParams, FindClusters, Network};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/clustering/edge")
}

fn tsv(p: &Path) -> Vec<BTreeMap<String, String>> {
    let text = std::fs::read_to_string(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    let mut lines = text.lines();
    let header: Vec<String> = lines.next().unwrap().split('\t').map(String::from).collect();
    lines.map(|l| header.iter().cloned().zip(l.split('\t').map(String::from)).collect()).collect()
}

fn sha(b: &[u8]) -> String {
    scio::sha256::hex(b)
}

fn lines(v: &[String]) -> Vec<u8> {
    v.iter().flat_map(|l| format!("{l}\n").into_bytes()).collect()
}

/// The files R's run of one setting wrote, recomputed from Rust's result: name -> bytes.
fn files_of(r: &FindClusters, barcodes: &[String], dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut out = BTreeMap::new();
    let _ = std::fs::remove_dir_all(dir);
    std::fs::create_dir_all(dir).unwrap();
    for name in cluster::output::write_labels(dir, barcodes, r).unwrap() {
        out.insert(name.clone(), std::fs::read(dir.join(&name)).unwrap());
    }
    out.insert("idents_levels.txt".into(), lines(&r.seurat_clusters.levels));
    let n = r.n_cells;
    let calls = r.resolutions.len();
    let col_major = |f: &dyn Fn(usize, usize) -> i32, rows: usize, cols: usize| -> Vec<i32> {
        let mut v = Vec::with_capacity(rows * cols);
        for c in 0..cols {
            for i in 0..rows {
                v.push(f(i, c));
            }
        }
        v
    };
    let opt = col_major(&|i, k| r.resolutions[k].optimizer.ids[i] as i32, n, calls);
    out.insert("stages/K7_optimizer_ids.npy".into(), scio::npy::encode_i32(&opt, &[n, calls]));
    let fin = col_major(&|i, k| r.resolutions[k].grouped.labels[i].map_or(-1, |v| v as i32), n, calls);
    out.insert("stages/K7_final_ids.npy".into(), scio::npy::encode_i32(&fin, &[n, calls]));
    // .Random.seed at each GroupSingletons call's entry and exit, after the script's
    // set.seed(20261114)
    let mut state = rcompat::RRng::new(20261114).random_seed();
    let mut rows: Vec<Vec<i32>> = Vec::new();
    for x in &r.resolutions {
        rows.push(state.clone());
        if let Some(s) = &x.grouped.rng_after {
            state = s.clone();
        }
        rows.push(state.clone());
    }
    let rng = col_major(&|i, j| rows[i][j], rows.len(), 626);
    out.insert("stages/rng.npy".into(), scio::npy::encode_i32(&rng, &[rows.len(), 626]));
    let mut t = String::from("row\tctx\tstage\tlabel\n");
    for k in 0..rows.len() {
        t.push_str(&format!("{k}\tchain\tK7\t{}\n", if k % 2 == 0 { "groupsingletons_entry" } else { "groupsingletons_exit" }));
    }
    out.insert("stages/rng.tsv".into(), t.into_bytes());
    out
}

/// The instrumented optimizer's files for one resolution, from Rust's optimizer.
fn optimizer_files(o: &cluster::OptimizerResult, res_name: &str, out: &mut BTreeMap<String, Vec<u8>>) {
    let n = o.n_nodes;
    let mut starts = String::from("start\titerations\tmodularity\n");
    for (s, st) in o.starts.iter().enumerate() {
        starts.push_str(&format!("{}\t{}\t{}\n", s + 1, st.iterations, scio::fmt17(st.modularity)));
    }
    out.insert(format!("optimizer/{res_name}/starts.tsv"), starts.into_bytes());
    let mut p = Vec::with_capacity(n * o.starts.len());
    for c in 0..n {
        for st in &o.starts {
            p.push(st.partition.as_ref().unwrap()[c] as i32);
        }
    }
    out.insert(format!("optimizer/{res_name}/partitions.npy"), scio::npy::encode_i32(&p, &[o.starts.len(), n]));
    let ids: Vec<String> = o.ids.iter().map(|v| v.to_string()).collect();
    out.insert(format!("optimizer/{res_name}/ids.txt"), lines(&ids));
}

#[test]
fn edge_cases_against_r() {
    let fx = fixture();
    let cases = tsv(&fx.join("cases.tsv"));
    let expected = tsv(&fx.join("expected.tsv"));
    let printed = tsv(&fx.join("printed.tsv"));
    assert!(cases.len() >= 50, "{} settings", cases.len());
    let tmp = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/test-tmp/cluster-edge");
    let (mut files_checked, mut errors_checked, mut pow2_reduced, mut ties_over_8, mut merged_order) = (0, 0, 0usize, 0usize, false);
    for c in &cases {
        let (case, setting) = (&c["case"], &c["setting"]);
        let input = cluster::input::read_snn_dir(&fx.join(case)).unwrap_or_else(|e| panic!("{case}: {e}"));
        let params = ClusterParams {
            modularity_fxn: c["modularity_fxn"].parse().unwrap(),
            resolutions: c["resolution"].split(',').map(|v| v.parse::<f64>().unwrap()).collect(),
            algorithm: c["algorithm"].parse().unwrap(),
            n_start: c["n_start"].parse().unwrap(),
            n_iter: c["n_iter"].parse().unwrap(),
            random_seed: c["random_seed"].parse().unwrap(),
            group_singletons: c["group_singletons"] == "TRUE",
            ..ClusterParams::default()
        };
        let label = format!("{case}/{setting}");
        let want: Vec<&BTreeMap<String, String>> = expected.iter().filter(|e| &e["case"] == case && &e["setting"] == setting).collect();
        let said = |kind: &str| -> Vec<String> { printed.iter().filter(|p| &p["case"] == case && &p["setting"] == setting && p["kind"] == kind).map(|p| p["text"].clone()).collect() };
        // the instrumented optimizer's files, from the optimizer alone (R ran it wherever R's own
        // optimizer finished, even where GroupSingletons then stopped)
        let mut got = BTreeMap::new();
        if let Some(net) = Network::from_snn(&input.snn, params.modularity_fxn) {
            for &res in &params.resolutions {
                let o = cluster::optimizer::run(&net, &params.optimizer(res), true).unwrap();
                optimizer_files(&o, &ClusterParams::res_name(res), &mut got);
                // reduced networks and subnetworks of 2^k nodes take nextInt's other branch: count
                // those draws where the input itself is not a power of two
                if !net.n_nodes.is_power_of_two() {
                    pow2_reduced += o.pow2_draws as usize;
                }
            }
        }
        let result = find_clusters(&input.snn, &params, 2, true);
        if !c["error"].is_empty() {
            let e = result.expect_err(&label).to_string();
            assert!(e.contains(&c["error"]), "{label}: R stopped with {:?}, Rust with {e:?}", c["error"]);
            assert_eq!(said("error"), vec![c["error"].clone()]);
            errors_checked += 1;
        } else {
            let r = result.unwrap_or_else(|e| panic!("{label}: {e}"));
            got.extend(files_of(&r, &input.barcodes, &tmp.join(case).join(setting)));
            // what FindClusters printed and its messages
            let rp: Vec<String> = r.resolutions.iter().flat_map(|x| {
                let mut p = x.printed();
                p.pop();
                p
            }).collect();
            assert_eq!(rp, said("printed"), "{label}: printed");
            let rm: Vec<String> = r.resolutions.iter().filter_map(|x| x.message()).collect();
            assert_eq!(rm, said("message"), "{label}: messages");
            let rs: Vec<String> = r.resolutions.iter().map(|x| x.grouped.singletons.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",")).collect();
            assert_eq!(rs, said("singletons"), "{label}: singletons");
            ties_over_8 += r.resolutions.iter().flat_map(|x| x.grouped.n_tied.iter()).filter(|&&t| t > 8).count();
            if case == "singleton_order" && setting == "modularity2" {
                // the second singleton joins the first one's cluster only because of that merge
                let g = &r.resolutions[0].grouped;
                merged_order = g.joined.len() == 2 && g.joined[0] == g.joined[1];
            }
        }
        for w in &want {
            let f = &w["file"];
            let bytes = got.get(f).unwrap_or_else(|| panic!("{label}: Rust has no {f}"));
            assert_eq!(sha(bytes), w["sha256"], "{label}: {f} ({} bytes in Rust's, {} in R's)", bytes.len(), w["bytes"]);
            files_checked += 1;
        }
        // nothing Rust wrote is missing from R's (for the files R writes)
        for f in got.keys() {
            assert!(want.iter().any(|w| &w["file"] == f) || f.ends_with(".tsv") && !f.contains('/') && f != "metadata.tsv",
                    "{label}: R wrote no {f}");
        }
    }
    assert!(files_checked > 650, "{files_checked} files");
    assert_eq!(errors_checked, 3);
    assert!(ties_over_8 >= 5, "{ties_over_8} singletons tied with more than 8 clusters");
    assert!(merged_order, "singleton_order: the second singleton should follow the first");
    assert!(pow2_reduced > 100, "{pow2_reduced} power-of-two draws on reduced networks");
}
