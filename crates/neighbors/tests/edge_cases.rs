//! Constructed edge cases (docs/PLAN.md C7), each against Seurat 5.5.1's FindNeighbors and
//! RcppAnnoy 0.0.23 run on the same input in the pinned R env (validation/R/neighbors_references.R
//! --mode edge, milestone M9), for all three metrics: duplicated cells (more identical cells than
//! k), margins of exactly 0, splits that fail and fall back to random sides, indexes no larger than
//! a leaf (n <= K), fewer cells than k.param, sizes that are powers of two, zero vectors, collinear
//! cells, and non-default n.trees, k.param and prune.SNN.
//!
//! The committed fixture holds each distinct input (`inputs/<sha256>.npy`), each case's parameters
//! and input (`cases.tsv`) and the size and SHA-256 of every file R wrote (`expected.tsv`): Annoy's index
//! file (`annoy.ann`), the items, the neighbour ids and distances, and both graphs' slots
//! (`dgCMatrix` p, i, x), all as `.npy` exactly as R's `npy_write` writes them. The R outputs
//! themselves are in data/neighbors/edge/ (listed in validation/manifest.toml). Matching digests
//! mean byte-identical files.
//!
//! The test also checks that the cases reach, with each metric, the branches they are meant for,
//! from the build's own counts: margins of exactly 0, splits retried, kept imbalanced and
//! randomised, and cells missing from their own neighbour lists.

use neighbors::{build_index, find_neighbors, items, search, Metric, NeighborParams};
use scio::npy::{encode_f64, encode_i32};
use scio::{Npy, Table};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/neighbors/edge")
}

fn sha(bytes: &[u8]) -> String {
    scio::sha256::hex(bytes)
}

/// Column-major copy of a row-major n x k array.
fn col_major<T: Copy>(rows: &[T], n: usize, k: usize) -> Vec<T> {
    (0..k).flat_map(|j| (0..n).map(move |i| (i, j))).map(|(i, j)| rows[i * k + j]).collect()
}

#[test]
fn edge_cases_match_r() {
    let cases = Table::read(&dir().join("cases.tsv")).unwrap();
    let expected = Table::read(&dir().join("expected.tsv")).unwrap();
    let col = |t: &Table, name: &str| t.col(name).unwrap_or_else(|| panic!("no column {name}"));
    let (ec, em, ef, eb, es) = (col(&expected, "case"), col(&expected, "metric"), col(&expected, "file"), col(&expected, "bytes"), col(&expected, "sha256"));
    let mut want: BTreeMap<(String, String), Vec<(String, usize, String)>> = BTreeMap::new();
    for r in &expected.rows {
        want.entry((r[ec].clone(), r[em].clone())).or_default().push((r[ef].clone(), r[eb].parse().unwrap(), r[es].clone()));
    }
    let mut compared = 0usize;
    let mut stats_seen: BTreeMap<(&str, &str), usize> = BTreeMap::new();
    let mut failures = Vec::new();
    for r in &cases.rows {
        let get = |name: &str| r[col(&cases, name)].clone();
        let case = get("case");
        let input = dir().join("inputs").join(format!("{}.npy", get("embeddings_sha256")));
        assert_eq!(sha(&std::fs::read(&input).unwrap()), get("embeddings_sha256"), "{case}: input differs from R's");
        let emb = Npy::read(&input).unwrap();
        let (n, f) = (emb.nrow(), emb.ncol());
        assert_eq!((n.to_string(), f.to_string()), (get("n"), get("f")));
        let data = emb.block_f64(0, n);
        for metric_name in ["euclidean", "cosine", "manhattan"] {
            let key = (case.clone(), metric_name.to_string());
            let files = want.get(&key).unwrap_or_else(|| panic!("{case} {metric_name}: no expected files"));
            let metric = Metric::parse(metric_name).unwrap();
            let params = NeighborParams {
                dims: (1..=f).collect(),
                k_param: get("k_param").parse().unwrap(),
                n_trees: get("n_trees").parse().unwrap(),
                metric,
                prune_snn: scio::parse17(&get("prune")).unwrap(),
                search_k: -1,
            };
            let mut got: BTreeMap<&str, Vec<u8>> = BTreeMap::new();
            // AnnoyBuildIndex and its file
            let it32 = items(&data, n, f, &params.dims);
            let index = build_index(&it32, n, f, metric, params.n_trees).unwrap();
            got.insert("stages/annoy.ann", index.to_bytes());
            let it: Vec<f64> = (0..n).flat_map(|i| index.get_item(i).iter().map(|&v| v as f64).collect::<Vec<_>>()).collect();
            got.insert("stages/K3_items.npy", encode_f64(&col_major(&it, n, f), &[n, f]));
            // AnnoySearch with FindNeighbors' k (the number of cells minus 1 when there are fewer)
            let k = if n < params.k_param { n - 1 } else { params.k_param };
            let (nn_idx, nn_dist, counts) = search(&index, &it32, k, -1, 2);
            assert!(counts.iter().all(|&c| c == k), "{case} {metric_name}: fewer than k neighbours");
            let ids: Vec<i32> = nn_idx.iter().map(|&i| i + 1).collect();
            got.insert("stages/K5_nn_index.npy", encode_i32(&col_major(&ids, n, k), &[n, k]));
            let d: Vec<f64> = nn_dist.iter().map(|&x| x as f64).collect();
            got.insert("stages/K5_nn_dist.npy", encode_f64(&col_major(&d, n, k), &[n, k]));
            let missing_self = (0..n).filter(|&c| !nn_idx[c * k..(c + 1) * k].contains(&(c as i32))).count();
            *stats_seen.entry((metric_name, "cells missing from their own list")).or_default() += missing_self;
            // the whole FindNeighbors call, with both graphs; R could not build them for one cell,
            // and the port refuses that case too
            let r_has_graphs = files.iter().any(|(fl, _, _)| fl == "graphs/snn_x.npy");
            match find_neighbors(&data, n, f, &params, 2) {
                Ok(res) => {
                    assert!(r_has_graphs, "{case} {metric_name}: R's FindNeighbors could not build the graphs, Rust did");
                    assert!(res.index.to_bytes() == got["stages/annoy.ann"] && res.nn_idx == nn_idx, "{case} {metric_name}: find_neighbors differs from its steps");
                    for (name, g) in [("nn", &res.nn), ("snn", &res.snn)] {
                        let p: Vec<i32> = g.p.iter().map(|&v| v as i32).collect();
                        let i: Vec<i32> = g.i.iter().map(|&v| v as i32).collect();
                        got.insert(if name == "nn" { "graphs/nn_p.npy" } else { "graphs/snn_p.npy" }, encode_i32(&p, &[p.len()]));
                        got.insert(if name == "nn" { "graphs/nn_i.npy" } else { "graphs/snn_i.npy" }, encode_i32(&i, &[i.len()]));
                        got.insert(if name == "nn" { "graphs/nn_x.npy" } else { "graphs/snn_x.npy" }, encode_f64(&g.x, &[g.x.len()]));
                    }
                }
                Err(e) => assert!(!r_has_graphs, "{case} {metric_name}: {e}"),
            }
            let st = index.stats();
            *stats_seen.entry((metric_name, "margins of exactly 0")).or_default() += st.zero_margins;
            *stats_seen.entry((metric_name, "splits accepted at the 2nd or 3rd attempt")).or_default() += st.accepted_at_attempt[1] + st.accepted_at_attempt[2];
            *stats_seen.entry((metric_name, "splits kept imbalanced")).or_default() += st.kept_imbalanced;
            *stats_seen.entry((metric_name, "splits with random sides")).or_default() += st.random_fallbacks;
            for (fl, bytes, digest) in files {
                let Some(b) = got.get(fl.as_str()) else {
                    failures.push(format!("{case} {metric_name} {fl}: not produced"));
                    continue;
                };
                compared += 1;
                if b.len() != *bytes || sha(b) != *digest {
                    failures.push(format!("{case} {metric_name} {fl}: {} bytes, R's {bytes}; digest {}", b.len(), if sha(b) == *digest { "equal" } else { "differs" }));
                }
            }
            for fl in got.keys() {
                if !files.iter().any(|(f2, _, _)| f2 == fl) {
                    failures.push(format!("{case} {metric_name} {fl}: R wrote none"));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{} of {compared} files differ:\n{}", failures.len(), failures.join("\n"));
    assert!(compared >= 700, "only {compared} files compared");
    assert_eq!(stats_seen.len(), 15);
    for ((metric, what), count) in &stats_seen {
        eprintln!("{metric}: {what}: {count}");
        assert!(*count > 0, "no case reached, with {metric}: {what}");
    }
}
