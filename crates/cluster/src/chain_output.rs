//! `scrust cluster`'s outputs from the chain (docs/PLAN.md C3, "Outputs"), written by an
//! [`Observer`] as each stage finishes, so the stages need not be kept:
//!
//! - `pca/`: RunPCA's result in `scrust pca`'s layout (`pca::output`);
//! - `neighbors/`: FindNeighbors' result in `scrust neighbors`' layout (`neighbors::output`);
//! - the labels' files (`seurat_clusters.tsv`, `clusters.tsv`, the levels, `cluster_ids.txt`,
//!   `nclusters.txt`, `optimizer_ids.tsv`; [`crate::output`]);
//! - `metadata.tsv`: `barcode`, `nCount_RNA` and `nFeature_RNA` (when known), each `res.`
//!   column, `seurat_clusters`, `silhouette_width` and `silhouette_other`, `neighborhood_purity`
//!   and `purity_maximum` (when computed), then the columns of `--metadata`;
//! - with `--keep-stages`: `sct/` in `scrust sctransform`'s layout (`sctransform::output`), and
//!   every random start's partition (`stages/<column>_partitions.npy`);
//! - with `--dump-stages <level>`: `stages/`, the stage dumps of M8's instrumented runs in their
//!   layout (context `chain`): S2a-S2d, S3, K3 (Annoy's items), K4 (`annoy.ann`), K5, K7, K8, K9
//!   (purity; its neighbour lists with `purity_lists`; the certificate's flagged pairs, which R
//!   has not), the `.Random.seed` checkpoints and the scalars, so `scrust validate cluster`
//!   compares them block by block.
//!
//! Each staged start writes what its stages make, and copies what it was given (`pca/` and
//! `neighbors/` from `scrust pca`'s and `scrust neighbors`' directories), so every route writes
//! the same files.

use crate::chain::Observer;
use crate::find_clusters::FindClusters;
use crate::metadata::RnaMetadata;
use crate::purity::Purity;
use crate::silhouette::Silhouette;
use neighbors::Neighbors;
use pca::Pca;
use scio::stages::{DumpLevel, StageRecorder};
use scio::{fmt17, Json};
use sctransform::Sct;
use std::io;
use std::path::{Path, PathBuf};

/// The context of the chain's stage dumps (M8's instrumented runs have one, `chain`).
pub const CTX: &str = "chain";

/// R's digest of a `.Random.seed` (`cl_rng_digest()` in `validation/R/lib/clustering.R`): the
/// SHA-256 of its integers joined by commas.
pub fn rng_digest(state: &[i32]) -> String {
    let text: Vec<String> = state.iter().map(|v| v.to_string()).collect();
    scio::sha256::hex(text.join(",").as_bytes())
}

/// Writes the chain's outputs into `out`.
pub struct Writer {
    pub out: PathBuf,
    pub barcodes: Vec<String>,
    /// The counts' features (id, name), when known: `sct/`'s genes and the loadings' rows.
    pub features: Option<Vec<(String, String)>>,
    pub rna: Option<RnaMetadata>,
    /// `--metadata`'s columns: names and, per column, the cells' values.
    pub extra: Option<(Vec<String>, Vec<Vec<String>>)>,
    pub keep_stages: bool,
    pub annoy: bool,
    /// Purity's neighbour lists in the stage dumps (R's `--purity-lists`).
    pub purity_lists: bool,
    pub threads: usize,
    /// FindNeighbors' parameters, for the stage dumps' scalars.
    pub neighbor_params: Option<neighbors::NeighborParams>,
    /// `seurat_clusters`' levels, once clustered.
    levels: Vec<String>,
    dumps: Option<StageRecorder>,
    /// What each stage contributes to `summary.json`'s `result` section.
    pub sct_result: Option<Json>,
    pub pca_result: Option<Json>,
    pub neighbors_result: Option<Json>,
    pub silhouette_result: Option<Json>,
    pub purity_result: Option<Json>,
    /// R's digests of `.Random.seed` after SCTransform and after RunPCA, and of `scale.data`.
    pub rng_after_sct: Option<String>,
    pub rng_after_pca: Option<String>,
    pub scale_data_sha256: Option<String>,
    /// The files written, relative to `out`.
    pub files: Vec<String>,
    /// Elapsed time (since the writer was made, [`Writer::started`]) and peak memory after each
    /// stage, for the run record.
    pub phases: Vec<(String, f64, Option<u64>)>,
    start: std::time::Instant,
}

/// Peak resident set size of this process so far (`VmHWM`), in KiB.
pub fn peak_rss_kib() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    s.lines().find_map(|l| l.strip_prefix("VmHWM:")).and_then(|v| v.trim().trim_end_matches("kB").trim().parse().ok())
}

fn mkdir(d: &Path) -> io::Result<()> {
    std::fs::create_dir_all(d).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", d.display())))
}

impl Writer {
    pub fn new(out: &Path, barcodes: Vec<String>, threads: usize) -> Writer {
        Writer {
            out: out.to_path_buf(),
            barcodes,
            features: None,
            rna: None,
            extra: None,
            keep_stages: false,
            annoy: false,
            purity_lists: false,
            threads,
            neighbor_params: None,
            levels: Vec::new(),
            dumps: None,
            sct_result: None,
            pca_result: None,
            neighbors_result: None,
            silhouette_result: None,
            purity_result: None,
            rng_after_sct: None,
            rng_after_pca: None,
            scale_data_sha256: None,
            files: Vec::new(),
            phases: Vec::new(),
            start: std::time::Instant::now(),
        }
    }

    /// When the writer was made: `phases` are timed from it.
    pub fn started(&self) -> std::time::Instant {
        self.start
    }

    /// Records stage dumps at `level` into `<out>/stages`.
    pub fn with_dumps(mut self, level: DumpLevel) -> Writer {
        let mut r = StageRecorder::new(&self.out.join("stages"), level);
        r.context(CTX);
        self.dumps = Some(r);
        self
    }

    pub fn dump_level(&self) -> Option<DumpLevel> {
        self.dumps.as_ref().map(|d| d.level)
    }

    /// Copies a stage directory's files into `<out>/<sub>` byte for byte (a staged start's
    /// `pca/` or `neighbors/`), with their `summary.json`'s `result` section for this run's.
    pub fn copy_stage(&mut self, from: &Path, sub: &str, files: &[&str]) -> io::Result<Option<Json>> {
        let d = self.out.join(sub);
        mkdir(&d)?;
        for f in files {
            let src = from.join(f);
            if src.is_file() {
                std::fs::copy(&src, d.join(f)).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", src.display())))?;
                self.files.push(format!("{sub}/{f}"));
            }
        }
        let summary = from.join("summary.json");
        Ok(std::fs::read_to_string(&summary).ok().and_then(|t| Json::parse(&t).ok()).and_then(|j| j.get("result").cloned()))
    }

    /// The labels' files and `metadata.tsv`, once the clustering (and the scores) are done.
    pub fn finish_labels(&mut self, r: &FindClusters, s: Option<&Silhouette>, p: Option<&Purity>) -> io::Result<()> {
        mkdir(&self.out)?;
        let written = crate::output::write_labels(&self.out, &self.barcodes, r)?;
        self.files.extend(written.into_iter().filter(|f| f != "metadata.tsv"));
        self.write_metadata(r, s, p)?;
        if self.keep_stages {
            let st = self.out.join("stages");
            mkdir(&st)?;
            self.files.extend(crate::output::write_partitions(&st, r)?.into_iter().map(|f| format!("stages/{f}")));
        }
        Ok(())
    }

    fn write_metadata(&mut self, r: &FindClusters, s: Option<&Silhouette>, p: Option<&Purity>) -> io::Result<()> {
        let n = self.barcodes.len();
        let mut names: Vec<String> = Vec::new();
        let mut cols: Vec<Vec<String>> = Vec::new();
        if let Some(rna) = &self.rna {
            let (a, b) = rna.texts();
            names.push("nCount_RNA".into());
            cols.push(a);
            names.push("nFeature_RNA".into());
            cols.push(b);
        }
        for x in &r.resolutions {
            names.push(x.column.name.clone());
            cols.push(x.column.labels().into_iter().map(String::from).collect());
        }
        names.push("seurat_clusters".into());
        cols.push(r.seurat_clusters.labels().into_iter().map(String::from).collect());
        if let Some(s) = s {
            names.push("silhouette_width".into());
            cols.push(s.width.iter().map(|&w| if w.is_nan() && !s.computed { "NA".to_string() } else { fmt17(w) }).collect());
            names.push("silhouette_other".into());
            cols.push(s.other_labels(&r.seurat_clusters));
        }
        if let Some(p) = p {
            names.push("neighborhood_purity".into());
            cols.push(p.purity.iter().map(|&v| if v.is_nan() && !p.computed { "NA".to_string() } else { fmt17(v) }).collect());
            names.push("purity_maximum".into());
            cols.push(p.maximum_labels(&r.seurat_clusters));
        }
        if let Some((en, ec)) = &self.extra {
            names.extend(en.iter().cloned());
            cols.extend(ec.iter().cloned());
        }
        let mut t = String::from("barcode");
        for h in &names {
            t.push('\t');
            t.push_str(h);
        }
        t.push('\n');
        for i in 0..n {
            t.push_str(&self.barcodes[i]);
            for c in &cols {
                t.push('\t');
                t.push_str(&c[i]);
            }
            t.push('\n');
        }
        let p = self.out.join("metadata.tsv");
        std::fs::write(&p, t).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", p.display())))?;
        self.files.push("metadata.tsv".into());
        Ok(())
    }

    /// Writes the stage dumps (`stages/index.tsv`, the arrays, `rng.npy`, `rng.tsv`,
    /// `stages.json`).
    pub fn finish_dumps(&mut self) -> io::Result<()> {
        if let Some(d) = self.dumps.take() {
            let names = d.writer.names();
            d.finish(&format!("scrust cluster {}", env!("CARGO_PKG_VERSION")))?;
            self.files.extend(names.into_iter().map(|f| format!("stages/{f}")));
            for f in ["index.tsv", "rng.npy", "rng.tsv", "stages.json"] {
                self.files.push(format!("stages/{f}"));
            }
        }
        Ok(())
    }

    /// `summary.json`'s `result` section: FindClusters' result (as `scrust cluster --snn` writes
    /// it), then each stage this run made or was given, R's digests and the scores.
    pub fn result_json(&self, r: &FindClusters) -> Json {
        let Json::Object(mut m) = crate::output::result_json(r) else { unreachable!() };
        let opt = |j: &Option<Json>| j.clone().unwrap_or(Json::Null);
        m.push(("sctransform".into(), opt(&self.sct_result)));
        m.push(("scale_data_sha256".into(), self.scale_data_sha256.clone().map(Json::Str).unwrap_or(Json::Null)));
        m.push(("pca".into(), opt(&self.pca_result)));
        m.push((
            "rng_after_sha256".into(),
            Json::obj(vec![
                ("SCTransform", self.rng_after_sct.clone().map(Json::Str).unwrap_or(Json::Null)),
                ("RunPCA", self.rng_after_pca.clone().map(Json::Str).unwrap_or(Json::Null)),
            ]),
        ));
        m.push(("neighbors".into(), opt(&self.neighbors_result)));
        m.push(("silhouette".into(), opt(&self.silhouette_result)));
        m.push(("purity".into(), opt(&self.purity_result)));
        Json::Object(m)
    }

    fn rec(&mut self) -> Option<&mut StageRecorder> {
        self.dumps.as_mut()
    }
}

impl Observer for Writer {
    fn sct(&mut self, sct: &Sct) -> io::Result<()> {
        let nf = self.features.as_ref().map_or(0, |f| f.len());
        self.sct_result = Some(sctransform::output::result_json(sct, nf));
        self.scale_data_sha256 = Some(sctransform::output::scale_data_sha256_r_order(&sct.residuals));
        self.rng_after_sct = Some(rng_digest(&sct.rng_exit));
        if self.keep_stages {
            let d = self.out.join("sct");
            mkdir(&d)?;
            let genes: Vec<sctransform::output::Gene> = self
                .features
                .as_ref()
                .expect("the counts' features")
                .iter()
                .map(|(id, name)| sctransform::output::Gene { id: id.clone(), name: name.clone() })
                .collect();
            let files = sctransform::output::write_outputs(&d, &self.barcodes, &genes, sct)?;
            self.files.extend(files.into_iter().map(|f| format!("sct/{f}")));
        }
        if let Some(r) = self.rec() {
            sctransform::record::record(r, CTX, &sct.setup, &sct.fits, sct.mean_cell_sum, &sct.reg, &sct.residuals, &sct.rng_entry, &sct.rng_exit)?;
        }
        Ok(())
    }

    fn pca(&mut self, input_rows: &[usize], input: &[f64], p: &Pca, entry: Option<&[i32]>, exit: &[i32]) -> io::Result<()> {
        let d = self.out.join("pca");
        mkdir(&d)?;
        let feats: Vec<(String, String)> = match &self.features {
            Some(f) => input_rows.iter().map(|&r| f[r].clone()).collect(),
            None => return Err(io::Error::other("the PCA's features are not known")),
        };
        pca::output::write_outputs(&d, &self.barcodes, &feats, p, exit)?;
        self.files.extend(pca::output::FILES.iter().map(|f| format!("pca/{f}")));
        self.pca_result = Some(pca::output::result_json(p, input_rows.len()));
        self.rng_after_pca = Some(rng_digest(exit));
        if let Some(r) = self.rec() {
            pca::record::record(r, CTX, input_rows, input, p, entry, exit)?;
        }
        Ok(())
    }

    fn neighbors(&mut self, r: &Neighbors, ncol: usize) -> io::Result<()> {
        let d = self.out.join("neighbors");
        mkdir(&d)?;
        let files = neighbors::output::write_outputs(&d, &self.barcodes, r, self.annoy, self.threads)?;
        self.files.extend(files.into_iter().map(|f| format!("neighbors/{f}")));
        self.neighbors_result = Some(neighbors::output::result_json(r, ncol));
        let stages = self.out.join("stages");
        let np = self.neighbor_params.clone();
        if let Some(rec) = self.rec() {
            let (n, k, f) = (r.n, r.k, r.index.f());
            // K3: the items Annoy holds (getItemsVector), floats as doubles, cells x dims
            let mut items = vec![0.0f64; n * f];
            for i in 0..n {
                let it = r.index.get_item(i);
                for z in 0..f {
                    items[i + z * n] = it[z] as f64;
                }
            }
            rec.writer.add_f64("K3_items", CTX, items, n, f, "K3");
            // K5: every cell's neighbours (1-based, as R holds them) and float distances
            let mut idx = vec![0i32; n * k];
            let mut dist = vec![0.0f64; n * k];
            for i in 0..n {
                for j in 0..k {
                    idx[i + j * n] = r.nn_idx[i * k + j] + 1;
                    dist[i + j * n] = r.nn_dist[i * k + j] as f64;
                }
            }
            rec.writer.add_i32("K5_nn_index", CTX, idx, n, k, "K5");
            rec.writer.add_f64("K5_nn_dist", CTX, dist, n, k, "K5");
            mkdir(&stages)?;
            r.index.save(&stages.join("annoy.ann"))?;
            let self_first = (0..n).filter(|&i| r.nn_idx[i * k] == i as i32).count();
            let self_missing = (0..n).filter(|&i| !r.nn_idx[i * k..(i + 1) * k].contains(&(i as i32))).count();
            let int = |v: usize| Json::Int(v as i64);
            if let Some(p) = &np {
                rec.scalar(CTX, "annoy_metric", Json::Str(p.metric.seurat_name().into()));
                rec.scalar(CTX, "annoy_n_trees_arg", Json::Int(p.n_trees as i64));
                rec.scalar(CTX, "knn_search_k", Json::Int(p.search_k as i64));
            }
            for (key, v) in [
                ("annoy_f", int(f)),
                ("annoy_n_items", int(n)),
                ("annoy_n_trees", int(r.index.n_trees())),
                ("annoy_file_bytes", int(r.index.n_nodes() * r.index.node_bytes())),
                ("annoy_node_bytes", int(r.index.node_bytes())),
                ("annoy_n_nodes", int(r.index.n_nodes())),
                ("knn_k", int(k)),
                ("knn_self_first", int(self_first)),
                ("knn_self_missing", int(self_missing)),
            ] {
                rec.scalar(CTX, key, v);
            }
        }
        if self.dumps.is_some() {
            self.files.push("stages/annoy.ann".into());
        }
        Ok(())
    }

    fn clusters(&mut self, fc: &FindClusters) -> io::Result<()> {
        self.levels = fc.seurat_clusters.levels.clone();
        if let Some(rec) = self.rec() {
            let n = fc.n_cells;
            let nres = fc.resolutions.len();
            let mut opt = Vec::with_capacity(n * nres);
            let mut fin = Vec::with_capacity(n * nres);
            for x in &fc.resolutions {
                opt.extend(x.optimizer.ids.iter().map(|&v| v as i32));
                fin.extend(x.grouped.labels.iter().map(|l| l.map_or(-1, |v| v as i32)));
            }
            rec.writer.add_i32("K7_optimizer_ids", CTX, opt, n, nres, "K7");
            rec.writer.add_i32("K7_final_ids", CTX, fin, n, nres, "K7");
            if let Some(s) = fc.rng_after() {
                rec.rng(CTX, "K7", "groupsingletons_exit", s.clone());
            }
            let p = &fc.params;
            for (key, v) in [
                ("optimizer_modularity_function", Json::Int(p.modularity_fxn as i64)),
                ("optimizer_resolution", Json::nums(&p.resolutions)),
                ("optimizer_algorithm", Json::Int(p.algorithm as i64)),
                ("optimizer_n_start", Json::Int(p.n_start as i64)),
                ("optimizer_n_iter", Json::Int(p.n_iter as i64)),
                ("optimizer_random_seed", Json::Int(p.random_seed as i64)),
            ] {
                rec.scalar(CTX, key, v);
            }
            let last = fc.resolutions.last().expect("a resolution");
            let o = &last.optimizer;
            let sing: Vec<String> = last.grouped.singletons.iter().map(|v| v.to_string()).collect();
            for (key, v) in [
                ("optimizer_communities", Json::Int(o.n_communities as i64)),
                ("singletons", Json::strs(&sing)),
                ("n_singletons", Json::Int(sing.len() as i64)),
            ] {
                rec.scalar(CTX, key, v);
            }
        }
        Ok(())
    }

    fn silhouette(&mut self, s: &Silhouette) -> io::Result<()> {
        let n = s.width.len();
        let self_levels = self.levels.clone();
        let computed = s.computed;
        self.silhouette_result = Some(Json::obj(vec![
            ("computed", Json::Bool(computed)),
            ("dims", Json::Vector(s.dims.iter().map(|&d| Json::Int(d as i64)).collect())),
            ("n_clusters", Json::Int(s.uclust.len() as i64)),
            ("cluster_var", Json::nums(&s.cluster_var)),
            ("widths_sha256", Json::Str(scio::sha256::hex(&s.width.iter().flat_map(|w| w.to_le_bytes()).collect::<Vec<u8>>()))),
        ]));
        if let Some(rec) = self.rec() {
            if computed {
                let nk = s.uclust.len();
                let nd = s.dims.len();
                rec.writer.add_f64("K8_silhouette_width", CTX, s.width.clone(), n, 1, "K8");
                rec.writer.add_i32("K8_silhouette_other", CTX, s.other.iter().map(|o| o.map_or(-1, |v| v as i32)).collect(), n, 1, "K8");
                // centroids: clusters x dims, column-major
                let mut c = vec![0.0; nk * nd];
                for k in 0..nk {
                    for d in 0..nd {
                        c[k + d * nk] = s.centroids[k * nd + d];
                    }
                }
                rec.writer.add_f64("K8_centroids", CTX, c, nk, nd, "K8");
                rec.writer.add_f64("K8_cluster_var", CTX, s.cluster_var.clone(), nk, 1, "K8");
            }
            rec.scalar(CTX, "silhouette_dims", Json::Int(s.dims.len() as i64));
            let lv: Vec<String> = s.uclust.iter().map(|&c| self_levels.get(c as usize).cloned().unwrap_or_default()).collect();
            rec.scalar(CTX, "silhouette_levels", Json::strs(&lv));
        }
        Ok(())
    }

    fn purity(&mut self, p: &Purity) -> io::Result<()> {
        let levels = self.levels.clone();
        self.purity_result = Some(purity_json(p, &levels, &self.barcodes));
        let lists = self.purity_lists;
        if let Some(rec) = self.rec() {
            if p.computed {
                let n = p.purity.len();
                let nk = p.uclust.len();
                rec.writer.add_f64("K9_dist50", CTX, p.dist_k.clone(), n, 1, "K9");
                rec.writer.add_i32("K9_neighbor_count", CTX, p.neighbor_count.iter().map(|&c| c as i32).collect(), n, 1, "K9");
                if lists {
                    if let Some((ptr, idx)) = &p.neighbor_lists {
                        rec.writer.add_i32("K9_neighbors_ptr", CTX, ptr.iter().map(|&v| v as i32).collect(), ptr.len(), 1, "K9");
                        rec.writer.add_i32("K9_neighbors_idx", CTX, idx.iter().map(|&v| v as i32 + 1).collect(), idx.len(), 1, "K9");
                    }
                }
                rec.writer.add_f64("K9_weights", CTX, p.weights.clone(), n, 1, "K9");
                // targets: cells x clusters, column-major
                let mut sums = vec![0.0; n * nk];
                for i in 0..n {
                    for c in 0..nk {
                        sums[i + c * n] = p.sums[i * nk + c];
                    }
                }
                rec.writer.add_f64("K9_sums", CTX, sums, n, nk, "K9");
                rec.writer.add_f64("K9_totals", CTX, p.totals.clone(), n, 1, "K9");
                rec.writer.add_f64("K9_purity", CTX, p.purity.clone(), n, 1, "K9");
                rec.writer.add_i32("K9_maximum", CTX, p.maximum.iter().map(|o| o.map_or(-1, |v| v as i32)).collect(), n, 1, "K9");
                // the certificate's flagged pairs (R has none): 1-based cells, the step (1:
                // findDistance, 3: findNeighbors), the kind (1: the median pair, on the radius by
                // construction; 2: other), the raw distance and the threshold
                let fl = &p.certificate.flags;
                if !fl.is_empty() {
                    let m = fl.len();
                    let mut pairs = vec![0i32; m * 2];
                    for (r, f) in fl.iter().enumerate() {
                        pairs[r] = f.i as i32 + 1;
                        pairs[r + m] = f.j as i32 + 1;
                    }
                    rec.writer.add_i32("K9_certificate_pairs", CTX, pairs, m, 2, "K9");
                    rec.writer.add_i32("K9_certificate_step", CTX, fl.iter().map(|f| if f.step == knn::purity::Step::KthDistance { 1 } else { 3 }).collect(), m, 1, "K9");
                    rec.writer.add_i32("K9_certificate_kind", CTX, fl.iter().map(|f| if f.kind == knn::purity::Kind::MedianPair { 1 } else { 2 }).collect(), m, 1, "K9");
                    rec.writer.add_f64("K9_certificate_raw", CTX, fl.iter().map(|f| f.raw).collect(), m, 1, "K9");
                    rec.writer.add_f64("K9_certificate_threshold", CTX, fl.iter().map(|f| f.threshold).collect(), m, 1, "K9");
                }
                rec.scalar(CTX, "purity_k", Json::Int(p.k_requested as i64));
                rec.scalar(CTX, "purity_radius", Json::Num(p.radius));
                rec.scalar(CTX, "purity_radius_is_median_of_dist50", Json::Bool(rcompat::quantile::median(&p.dist_k).to_bits() == p.radius.to_bits()));
                let lv: Vec<String> = p.uclust.iter().map(|&c| levels.get(c as usize).cloned().unwrap_or_default()).collect();
                rec.scalar(CTX, "purity_levels", Json::strs(&lv));
                rec.scalar(CTX, "purity_neighbor_lists", Json::Bool(lists && p.neighbor_lists.is_some()));
                rec.scalar(CTX, "purity_mean_neighbors", Json::Num(p.mean_neighbors()));
                rec.scalar(CTX, "purity_search", Json::Str("exhaustive, with knncolle's arithmetic".into()));
                rec.scalar(CTX, "purity_k_used", Json::Int(p.k as i64));
                rec.scalar(CTX, "purity_certificate_flagged", Json::Int(fl.len() as i64));
                rec.scalar(CTX, "purity_certificate_flagged_median_pair", Json::Int(p.certificate.n_flagged_kind(knn::purity::Kind::MedianPair) as i64));
                rec.scalar(CTX, "purity_certificate_flagged_other", Json::Int(p.certificate.n_flagged_kind(knn::purity::Kind::Other) as i64));
            }
        }
        Ok(())
    }

    fn phase(&mut self, name: &str) {
        self.phases.push((name.to_string(), self.start.elapsed().as_secs_f64(), peak_rss_kib()));
    }

}

/// `summary.json`'s `result.purity`: what was computed, its summary, and the certificate (every
/// flagged pair, with how close it lies to its threshold).
pub fn purity_json(p: &Purity, levels: &[String], barcodes: &[String]) -> Json {
    let num = |v: f64| if v.is_finite() { Json::Num(v) } else { Json::Null };
    let opt = |v: Option<f64>| v.map(num).unwrap_or(Json::Null);
    let c = &p.certificate;
    let mut items = vec![
        ("computed", Json::Bool(p.computed)),
        ("dims", Json::Vector(p.dims.iter().map(|&d| Json::Int(d as i64)).collect())),
        ("k", Json::Int(p.k_requested as i64)),
        ("k_used", Json::Int(p.k as i64)),
        ("weighted", Json::Bool(p.weighted)),
        ("n_clusters", Json::Int(p.uclust.len() as i64)),
        ("levels", Json::strs(&p.uclust.iter().map(|&c| levels.get(c as usize).cloned().unwrap_or_default()).collect::<Vec<_>>())),
    ];
    if p.computed {
        let flags: Vec<Json> = c
            .flags
            .iter()
            .map(|f| {
                let bc = |i: u32| Json::Str(barcodes.get(i as usize).cloned().unwrap_or_default());
                Json::obj(vec![
                    ("step", Json::Str(f.step.name().into())),
                    ("kind", Json::Str(f.kind.name().into())),
                    ("cell", Json::Int(f.i as i64 + 1)),
                    ("other", Json::Int(f.j as i64 + 1)),
                    ("cell_barcode", bc(f.i)),
                    ("other_barcode", bc(f.j)),
                    ("raw", Json::Num(f.raw)),
                    ("threshold", Json::Num(f.threshold)),
                    ("ulps_from_threshold", Json::Num(f.ulps())),
                    ("band_end", Json::Num(f.band)),
                    ("distance_minus_radius", if f.step == knn::purity::Step::Radius { Json::Num(f.raw.sqrt() - p.radius) } else { Json::Null }),
                ])
            })
            .collect();
        let mg = p.margin.as_ref();
        let t2 = p.radius * p.radius;
        items.extend([
            ("radius", Json::Num(p.radius)),
            ("radius_raw", Json::Num(t2)),
            ("mean_neighbors", Json::Num(p.mean_neighbors())),
            ("search", Json::Str("exhaustive, with knncolle's arithmetic (the result is KMKNN's unless the certificate flags a pair)".into())),
            ("purities_sha256", Json::Str(scio::sha256::hex(&p.purity.iter().flat_map(|w| w.to_le_bytes()).collect::<Vec<u8>>()))),
            (
                "certificate",
                Json::obj(vec![
                    ("flagged_pairs", Json::Int(c.n_flagged() as i64)),
                    ("flagged_kth_distance", Json::Int(c.n_flagged_step(knn::purity::Step::KthDistance) as i64)),
                    ("flagged_radius", Json::Int(c.n_flagged_step(knn::purity::Step::Radius) as i64)),
                    ("flagged_median_pair", Json::Int(c.n_flagged_kind(knn::purity::Kind::MedianPair) as i64)),
                    ("flagged_other", Json::Int(c.n_flagged_kind(knn::purity::Kind::Other) as i64)),
                    (
                        "treatment",
                        Json::Str("reported, not failing: the exhaustive result is kept (M12's decision); a median pair lies on the radius by construction, and R's KMKNN kept every one M12 measured".into()),
                    ),
                    ("flagged_cells", Json::Int({
                        let mut v: Vec<u32> = c.flags.iter().flat_map(|f| [f.i, f.j]).collect();
                        v.sort_unstable();
                        v.dedup();
                        v.len() as i64
                    })),
                    ("pairs", Json::List(flags)),
                    ("dims", mg.map(|m| Json::Int(m.dims as i64)).unwrap_or(Json::Null)),
                    ("gamma", mg.map(|m| Json::Num(m.gamma)).unwrap_or(Json::Null)),
                    ("eps", mg.map(|m| Json::Num(m.eps)).unwrap_or(Json::Null)),
                    ("kappa", mg.map(|m| Json::Num(m.kappa)).unwrap_or(Json::Null)),
                    ("bound", mg.map(|m| Json::Num(m.bound)).unwrap_or(Json::Null)),
                    ("radius_band_raw", num(c.radius_band)),
                    ("radius_band_ulps", num(c.radius_band / knn::purity::ulp(t2.max(f64::MIN_POSITIVE)))),
                    ("radius_gap_below_raw", opt(c.radius_gap_below)),
                    ("radius_gap_below_ulps", opt(c.radius_gap_below.map(|g| g / knn::purity::ulp(t2.max(f64::MIN_POSITIVE))))),
                    ("radius_gap_above_raw", opt(c.radius_gap_above)),
                    ("kth_band_max_raw", num(c.kth_band_max)),
                    ("kth_band_min_raw", num(c.kth_band_min)),
                    ("kth_gap_ratio", opt(c.kth_gap_ratio)),
                ]),
            ),
        ]);
    }
    Json::obj(items)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_s_rng_digest() {
        // digest::digest("1,2,3", algo = "sha256", serialize = FALSE)
        assert_eq!(rng_digest(&[1, 2, 3]), scio::sha256::hex(b"1,2,3"));
    }
}
