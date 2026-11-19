//! The output manifest (docs/PLAN.md C3f, milestone M13): the `files` section of `summary.json`,
//! in one schema for every scrust command that writes `summary.json`. It lists every file of the
//! output directory with its role, its scale, its format and shape, where its row and column
//! identifiers are, its size and its SHA-256, and records how the output was made (provenance).
//! The schema is versioned ([`SCHEMA`], [`VERSION`]) and documented in
//! `docs/output_manifest.md`, which the loaders (C3g, M14) are written against.
//!
//! Two rules keep normalised data from being taken for raw counts ([`naming_problems`]):
//! - a file named as a 10x reader looks for (`matrix.mtx`, `barcodes.tsv`, `features.tsv`,
//!   `genes.tsv`, gzipped or not) holds raw counts or their cells or features, and shares its
//!   directory only with such files: so no directory holding normalised data has a 10x name;
//! - no other file's name ends in `matrix.mtx` or `matrix.mtx.gz` (scanpy's `read_10x_mtx` with
//!   `prefix=` reads `<prefix>matrix.mtx.gz`).

use crate::json::Json;
use std::collections::BTreeSet;
use std::io;
use std::path::Path;

/// The schema's name and version, written in every manifest.
pub const SCHEMA: &str = "scrust-files";
pub const VERSION: i64 = 1;
/// Where the schema is described, relative to the repository.
pub const DOCUMENTATION: &str = "docs/output_manifest.md";

/// The file names 10x readers look for (Seurat's `Read10X`, scanpy's `read_10x_mtx`).
pub const TENX_NAMES: [&str; 8] = ["matrix.mtx", "matrix.mtx.gz", "barcodes.tsv", "barcodes.tsv.gz", "features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"];

/// The roles of a 10x directory of raw counts: the matrix, its cells, its features.
pub const RAW_ROLES: [&str; 3] = ["raw_counts", "raw_counts_cells", "raw_counts_features"];

/// The scales a numeric file can have. `raw_counts` is the only one a counts reader may take;
/// `none` marks files that hold no matrix of values (tables of mixed columns, labels, text).
pub const SCALES: [&str; 12] = [
    "raw_counts",
    "corrected_counts",
    "log1p_corrected_counts",
    "residuals",
    "pca_embeddings",
    "pca_loadings",
    "pca_stdev",
    "graph_nn",
    "graph_snn",
    "neighbour_indices",
    "neighbour_distances",
    "none",
];

/// The formats a file can have.
pub const FORMATS: [&str; 6] = ["matrix_market", "tsv", "npy", "json", "text", "annoy"];

/// The role of a nested output's `summary.json` (M15): a subdirectory that is itself a scrust
/// output, described by its own manifest. Its files are that manifest's, not the parent's: the
/// parent lists only the nested `summary.json` (with its size and SHA-256), and the checks recurse.
pub const NESTED_OUTPUT: &str = "nested_output";

/// How deep nested outputs may go before [`Manifest::check`] stops recursing (scrust writes one
/// level: `scrust doubletfinder`'s `clusters/` and `singlets_clusters/`).
pub const MAX_NESTING: usize = 8;

/// Whether a listed path stays inside the output it is listed in: relative, with no `..`, `.` or
/// empty component. No listed file, and above all no nested output, may escape its output.
pub fn path_inside(path: &str) -> bool {
    !path.starts_with('/') && path.split('/').all(|c| !c.is_empty() && c != "." && c != "..")
}

/// A column of the file holding identifiers: its header name, or its 1-based number when the file
/// has no header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Column {
    Name(String),
    Number(usize),
}

/// Where the identifiers of a file's rows or columns are.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ids {
    /// What they identify: `cells`, `features`, `variable_features`, `components`, `clusters`,
    /// `neighbours`, `starts`, `resolutions`.
    pub kind: String,
    /// The file that lists them, relative to the output directory (it may be the entry itself).
    pub file: String,
    pub column: Column,
}

impl Ids {
    /// Identifiers in the column `name` of a file with a header.
    pub fn named(kind: &str, file: &str, name: &str) -> Ids {
        Ids { kind: kind.into(), file: file.into(), column: Column::Name(name.into()) }
    }

    /// Identifiers in column `n` (1-based) of a file without a header (one per line when `n` is 1).
    pub fn numbered(kind: &str, file: &str, n: usize) -> Ids {
        Ids { kind: kind.into(), file: file.into(), column: Column::Number(n) }
    }

    fn to_json(&self) -> Json {
        Json::obj(vec![
            ("kind", Json::Str(self.kind.clone())),
            ("file", Json::Str(self.file.clone())),
            (
                "column",
                match &self.column {
                    Column::Name(n) => Json::Str(n.clone()),
                    Column::Number(k) => Json::Int(*k as i64),
                },
            ),
        ])
    }

    fn from_json(j: &Json) -> Result<Ids, String> {
        let kind = j.get("kind").and_then(Json::as_str).ok_or("identifiers without a kind")?.to_string();
        let file = j.get("file").and_then(Json::as_str).ok_or("identifiers without a file")?.to_string();
        let column = match j.get("column") {
            Some(Json::Str(s)) => Column::Name(s.clone()),
            Some(Json::Int(k)) if *k >= 1 => Column::Number(*k as usize),
            _ => return Err("identifiers without a column".into()),
        };
        Ok(Ids { kind, file, column })
    }
}

/// One file of an output directory.
#[derive(Debug, Clone, PartialEq)]
pub struct Entry {
    /// Relative to the output directory, with `/` between directories.
    pub path: String,
    /// What the file is (docs/output_manifest.md lists the roles).
    pub role: String,
    /// The scale of its values ([`SCALES`]).
    pub scale: String,
    /// [`FORMATS`]; `gzip` when the file is gzipped; `header` for tables (`tsv`): whether the
    /// first line names the columns.
    pub format: String,
    pub gzip: bool,
    pub header: Option<bool>,
    /// Matrix Market and tables: rows and columns (a table's data rows, without its header);
    /// `.npy`: the array's shape. `None` for other formats.
    pub shape: Option<Vec<usize>>,
    /// Matrix Market: the stored entries.
    pub nnz: Option<usize>,
    /// Where the rows' and the columns' identifiers are; `None` when there are none, or (columns
    /// of a table) when the file's own header names them.
    pub rows: Option<Ids>,
    pub columns: Option<Ids>,
    pub description: String,
    /// Size and SHA-256 of the file; `None` only for `summary.json` itself.
    pub bytes: Option<u64>,
    pub sha256: Option<String>,
}

impl Entry {
    pub fn new(path: &str, role: &str, scale: &str, format: &str, description: &str) -> Entry {
        Entry {
            path: path.into(),
            role: role.into(),
            scale: scale.into(),
            format: format.into(),
            gzip: path.ends_with(".gz"),
            header: None,
            shape: None,
            nnz: None,
            rows: None,
            columns: None,
            description: description.into(),
            bytes: None,
            sha256: None,
        }
    }

    pub fn header(mut self, h: bool) -> Entry {
        self.header = Some(h);
        self
    }

    pub fn shape(mut self, s: &[usize]) -> Entry {
        self.shape = Some(s.to_vec());
        self
    }

    pub fn nnz(mut self, n: usize) -> Entry {
        self.nnz = Some(n);
        self
    }

    pub fn rows(mut self, ids: Ids) -> Entry {
        self.rows = Some(ids);
        self
    }

    pub fn columns(mut self, ids: Ids) -> Entry {
        self.columns = Some(ids);
        self
    }

    /// The same entry under the directory `dir` (its identifier files too).
    pub fn under(mut self, dir: &str) -> Entry {
        if dir.is_empty() {
            return self;
        }
        let pre = |p: &str| format!("{}/{p}", dir.trim_end_matches('/'));
        self.path = pre(&self.path);
        for ids in [&mut self.rows, &mut self.columns].into_iter().flatten() {
            ids.file = pre(&ids.file);
        }
        self
    }

    /// The file's name, without its directories.
    pub fn file_name(&self) -> &str {
        self.path.rsplit('/').next().unwrap_or(&self.path)
    }

    /// The directory holding it (`""` at the top).
    pub fn dir(&self) -> &str {
        match self.path.rfind('/') {
            Some(k) => &self.path[..k],
            None => "",
        }
    }

    pub fn is_raw(&self) -> bool {
        RAW_ROLES.contains(&self.role.as_str())
    }

    pub fn to_json(&self) -> Json {
        let opt_ids = |i: &Option<Ids>| i.as_ref().map(Ids::to_json).unwrap_or(Json::Null);
        Json::obj(vec![
            ("path", Json::Str(self.path.clone())),
            ("role", Json::Str(self.role.clone())),
            ("scale", Json::Str(self.scale.clone())),
            ("format", Json::Str(self.format.clone())),
            ("gzip", Json::Bool(self.gzip)),
            ("header", self.header.map(Json::Bool).unwrap_or(Json::Null)),
            ("shape", self.shape.as_ref().map(|s| Json::Vector(s.iter().map(|&v| Json::Int(v as i64)).collect())).unwrap_or(Json::Null)),
            ("nnz", self.nnz.map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
            ("rows", opt_ids(&self.rows)),
            ("columns", opt_ids(&self.columns)),
            ("description", Json::Str(self.description.clone())),
            ("bytes", self.bytes.map(|v| Json::Int(v as i64)).unwrap_or(Json::Null)),
            ("sha256", self.sha256.clone().map(Json::Str).unwrap_or(Json::Null)),
        ])
    }

    pub fn from_json(j: &Json) -> Result<Entry, String> {
        let s = |k: &str| j.get(k).and_then(Json::as_str).map(String::from).ok_or_else(|| format!("an entry without {k}"));
        let path = s("path")?;
        let ids = |k: &str| match j.get(k) {
            None | Some(Json::Null) => Ok(None),
            Some(x) => Ids::from_json(x).map(Some).map_err(|e| format!("{path}: {k}: {e}")),
        };
        let shape = match j.get("shape") {
            Some(Json::Vector(v) | Json::List(v)) => Some(v.iter().map(|x| x.as_f64().map(|f| f as usize).ok_or_else(|| format!("{path}: bad shape"))).collect::<Result<Vec<_>, _>>()?),
            _ => None,
        };
        Ok(Entry {
            role: s("role")?,
            scale: s("scale")?,
            format: s("format")?,
            gzip: matches!(j.get("gzip"), Some(Json::Bool(true))),
            header: match j.get("header") {
                Some(Json::Bool(b)) => Some(*b),
                _ => None,
            },
            shape,
            nnz: j.get("nnz").and_then(Json::as_f64).map(|v| v as usize),
            rows: ids("rows")?,
            columns: ids("columns")?,
            description: s("description").unwrap_or_default(),
            bytes: j.get("bytes").and_then(Json::as_f64).map(|v| v as u64),
            sha256: j.get("sha256").and_then(Json::as_str).map(String::from),
            path,
        })
    }
}

/// A file a command read, with its size and SHA-256.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputFile {
    pub name: String,
    pub bytes: u64,
    pub sha256: String,
}

/// One input of a command: the option that named it, its path, and the files read from it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Input {
    pub option: String,
    pub path: String,
    pub files: Vec<InputFile>,
}

impl Input {
    /// `files` (paths) digested, named relative to `path` when under it.
    pub fn digest(option: &str, path: &Path, files: &[&Path]) -> io::Result<Input> {
        let mut out = Vec::with_capacity(files.len());
        for f in files {
            let (bytes, sha256) = crate::sha256::file(f)?;
            let name = f.strip_prefix(path).unwrap_or(f).display().to_string();
            out.push(InputFile { name, bytes, sha256 });
        }
        Ok(Input { option: option.into(), path: path.display().to_string(), files: out })
    }

    fn to_json(&self) -> Json {
        Json::obj(vec![
            ("option", Json::Str(self.option.clone())),
            ("path", Json::Str(self.path.clone())),
            (
                "files",
                Json::List(
                    self.files
                        .iter()
                        .map(|f| Json::obj(vec![("name", Json::Str(f.name.clone())), ("bytes", Json::Int(f.bytes as i64)), ("sha256", Json::Str(f.sha256.clone()))]))
                        .collect(),
                ),
            ),
        ])
    }

    fn from_json(j: &Json) -> Result<Input, String> {
        let files = match j.get("files") {
            Some(Json::List(v) | Json::Vector(v)) => v
                .iter()
                .map(|f| {
                    Ok(InputFile {
                        name: f.get("name").and_then(Json::as_str).ok_or("input file without a name")?.into(),
                        bytes: f.get("bytes").and_then(Json::as_f64).ok_or("input file without bytes")? as u64,
                        sha256: f.get("sha256").and_then(Json::as_str).ok_or("input file without sha256")?.into(),
                    })
                })
                .collect::<Result<Vec<_>, String>>()?,
            _ => Vec::new(),
        };
        Ok(Input {
            option: j.get("option").and_then(Json::as_str).unwrap_or("").into(),
            path: j.get("path").and_then(Json::as_str).unwrap_or("").into(),
            files,
        })
    }
}

/// How the output was made.
#[derive(Debug, Clone, PartialEq)]
pub struct Provenance {
    pub scrust_version: String,
    /// The command line as given (`scrust` and its arguments).
    pub command_line: Vec<String>,
    /// The command as the normalised matrices' comment lines give it: the command line without the
    /// options that change no value of any file (`--out`, `--threads`, ...), so that the files do
    /// not depend on them.
    pub command: String,
    pub inputs: Vec<Input>,
    /// The original software the outputs reproduce (R functions and packages, versions).
    pub reproduces: Vec<String>,
    pub parameters: Json,
    pub seeds: Json,
}

impl Provenance {
    fn to_json(&self) -> Json {
        Json::obj(vec![
            ("scrust_version", Json::Str(self.scrust_version.clone())),
            ("command_line", Json::strs(&self.command_line)),
            ("command", Json::Str(self.command.clone())),
            ("inputs", Json::List(self.inputs.iter().map(Input::to_json).collect())),
            ("reproduces", Json::strs(&self.reproduces)),
            ("parameters", self.parameters.clone()),
            ("seeds", self.seeds.clone()),
        ])
    }
}

/// The manifest of one output directory.
#[derive(Debug, Clone, PartialEq)]
pub struct Manifest {
    /// The scrust command: `cluster`, `sctransform`, `pca`, `neighbors`, `doubletfinder`.
    pub command: String,
    pub provenance: Provenance,
    pub entries: Vec<Entry>,
}

/// Every file under `dir`, relative to it with `/` separators, sorted.
pub fn list_files(dir: &Path) -> io::Result<Vec<String>> {
    fn walk(base: &Path, d: &Path, out: &mut Vec<String>) -> io::Result<()> {
        let rd = std::fs::read_dir(d).map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", d.display())))?;
        for e in rd {
            let e = e?;
            let p = e.path();
            let ft = std::fs::metadata(&p)?.file_type();
            if ft.is_dir() {
                walk(base, &p, out)?;
            } else {
                let rel = p.strip_prefix(base).unwrap_or(&p);
                out.push(rel.components().map(|c| c.as_os_str().to_string_lossy().into_owned()).collect::<Vec<_>>().join("/"));
            }
        }
        Ok(())
    }
    let mut out = Vec::new();
    walk(dir, dir, &mut out)?;
    out.sort();
    Ok(out)
}

/// The naming rules' violations among `entries` (none means the rules hold): see the module's
/// documentation.
pub fn naming_problems(entries: &[Entry]) -> Vec<String> {
    let mut out = Vec::new();
    let tenx_dirs: BTreeSet<&str> = entries.iter().filter(|e| TENX_NAMES.contains(&e.file_name())).map(|e| e.dir()).collect();
    for e in entries {
        let name = e.file_name();
        if TENX_NAMES.contains(&name) && !e.is_raw() {
            out.push(format!("{}: a 10x reader's file name on a file that is not raw counts (role {}, scale {})", e.path, e.role, e.scale));
        }
        if tenx_dirs.contains(e.dir()) && !e.is_raw() {
            out.push(format!(
                "{}: in {}, a directory with 10x file names, but not raw counts (role {}, scale {})",
                e.path,
                if e.dir().is_empty() { "the top directory" } else { e.dir() },
                e.role,
                e.scale
            ));
        }
        if !e.is_raw() && (name.ends_with("matrix.mtx") || name.ends_with("matrix.mtx.gz")) {
            out.push(format!("{}: a name ending in matrix.mtx(.gz) on a file that is not raw counts", e.path));
        }
        if e.is_raw() && e.role == "raw_counts" && e.scale != "raw_counts" {
            out.push(format!("{}: role raw_counts with scale {}", e.path, e.scale));
        }
    }
    out
}

impl Manifest {
    pub fn new(command: &str, provenance: Provenance) -> Manifest {
        Manifest { command: command.into(), provenance, entries: Vec::new() }
    }

    pub fn push(&mut self, e: Entry) {
        self.entries.push(e);
    }

    /// Fills the shapes the entries lack, from the files in `dir`: a table's data rows and columns
    /// (gzipped or not), an `.npy` array's shape from its header, a Matrix Market file's size line
    /// (decompressing it: callers give the shapes of large matrices they wrote).
    pub fn fill_shapes(&mut self, dir: &Path) -> io::Result<()> {
        for e in self.entries.iter_mut().filter(|e| e.shape.is_none()) {
            let p = dir.join(&e.path);
            match e.format.as_str() {
                "tsv" => {
                    let bytes = crate::tenx::read_maybe_gz(&p)?;
                    let text = String::from_utf8_lossy(&bytes);
                    let mut lines = text.lines();
                    let first = lines.next();
                    let ncol = first.map_or(0, |l| l.split('\t').count());
                    let nrow = lines.count() + usize::from(first.is_some() && e.header != Some(true));
                    e.shape = Some(vec![nrow, ncol]);
                }
                "npy" => e.shape = Some(crate::npy::Npy::read_header(&p).map_err(|x| io::Error::other(x.to_string()))?.1),
                "matrix_market" => {
                    let entries = crate::mtx::read_entries(&p)?;
                    e.shape = Some(vec![entries.nrow, entries.ncol]);
                    e.nnz = Some(entries.vals.len());
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// The entry for `summary.json`, the manifest's own file (it has no size or checksum).
    pub fn summary_entry() -> Entry {
        Entry::new("summary.json", "summary", "none", "json", "the run's summary: its results (result), how it ran (run) and this manifest (files)")
    }

    /// Fills every entry's size and SHA-256 from the files in `dir` (but `summary.json`'s), on up
    /// to `threads` threads; the result does not depend on `threads`. An entry that has both
    /// already keeps them: the command that wrote the file digested it as it wrote it.
    pub fn fill_digests(&mut self, dir: &Path, threads: usize) -> io::Result<()> {
        let todo: Vec<usize> = (0..self.entries.len()).filter(|&k| self.entries[k].path != "summary.json" && !(self.entries[k].bytes.is_some() && self.entries[k].sha256.is_some())).collect();
        let paths: Vec<std::path::PathBuf> = todo.iter().map(|&k| dir.join(&self.entries[k].path)).collect();
        let results: Vec<io::Result<(u64, String)>> = if threads <= 1 || paths.len() <= 1 {
            paths.iter().map(|p| crate::sha256::file(p)).collect()
        } else {
            // the largest files first, each on whichever thread is free; results placed by index
            let mut order: Vec<usize> = (0..paths.len()).collect();
            let size = |p: &Path| std::fs::metadata(p).map(|m| m.len()).unwrap_or(0);
            order.sort_by_key(|&k| std::cmp::Reverse(size(&paths[k])));
            let next = std::sync::atomic::AtomicUsize::new(0);
            let slots: Vec<std::sync::Mutex<Option<io::Result<(u64, String)>>>> = paths.iter().map(|_| std::sync::Mutex::new(None)).collect();
            std::thread::scope(|s| {
                for _ in 0..threads.min(paths.len()) {
                    s.spawn(|| loop {
                        let t = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if t >= order.len() {
                            break;
                        }
                        let k = order[t];
                        *slots[k].lock().unwrap() = Some(crate::sha256::file(&paths[k]));
                    });
                }
            });
            slots.into_iter().map(|s| s.into_inner().unwrap().expect("every file digested")).collect()
        };
        for (k, r) in todo.into_iter().zip(results) {
            let (bytes, sha256) = r?;
            self.entries[k].bytes = Some(bytes);
            self.entries[k].sha256 = Some(sha256);
        }
        Ok(())
    }

    /// Checks the manifest against the directory it describes: every file there is listed, every
    /// listed file exists with its recorded size and SHA-256 (`verify_digests`), each path is
    /// listed once, every role, scale and format is known, every identifier file is listed, and
    /// the naming rules hold; every listed path stays inside `dir`. A nested output is checked
    /// against its own manifest, recursively, its problems reported under its directory. Returns
    /// the problems found.
    pub fn check(&self, dir: &Path, verify_digests: bool) -> io::Result<Vec<String>> {
        self.check_at(dir, verify_digests, 0)
    }

    fn check_at(&self, dir: &Path, verify_digests: bool, depth: usize) -> io::Result<Vec<String>> {
        let mut out = Vec::new();
        let listed: BTreeSet<&str> = self.entries.iter().map(|e| e.path.as_str()).collect();
        if listed.len() != self.entries.len() {
            out.push("a path is listed more than once".into());
        }
        // a nested output's files are its own manifest's (checked below)
        let nested: Vec<String> = self.entries.iter().filter(|e| e.role == NESTED_OUTPUT && path_inside(&e.path)).map(|e| format!("{}/", e.dir())).collect();
        for f in list_files(dir)? {
            if !listed.contains(f.as_str()) && !nested.iter().any(|d| d != "/" && f.starts_with(d.as_str())) {
                out.push(format!("{f}: in the directory, not in the manifest"));
            }
        }
        for e in &self.entries {
            if !path_inside(&e.path) {
                out.push(format!("{}: not a relative path inside the output (absolute, or with a `..`, `.` or empty component)", e.path));
                continue;
            }
            let p = dir.join(&e.path);
            if e.path == "summary.json" {
                // the manifest's own file: no size or checksum (it holds them)
                continue;
            }
            if !p.is_file() {
                out.push(format!("{}: listed, not in the directory", e.path));
                continue;
            }
            if !SCALES.contains(&e.scale.as_str()) {
                out.push(format!("{}: unknown scale {}", e.path, e.scale));
            }
            if !FORMATS.contains(&e.format.as_str()) {
                out.push(format!("{}: unknown format {}", e.path, e.format));
            }
            if e.role.is_empty() {
                out.push(format!("{}: no role", e.path));
            }
            for ids in [&e.rows, &e.columns].into_iter().flatten() {
                if !listed.contains(ids.file.as_str()) {
                    out.push(format!("{}: its identifiers' file {} is not listed", e.path, ids.file));
                }
            }
            match (e.bytes, &e.sha256) {
                (Some(b), Some(s)) if verify_digests => {
                    let (b2, s2) = crate::sha256::file(&p)?;
                    if b != b2 || *s != s2 {
                        out.push(format!("{}: {b2} bytes, sha256 {s2}; the manifest says {b} bytes, {s}", e.path));
                    }
                }
                (Some(_), Some(_)) => {}
                _ => out.push(format!("{}: no size or checksum", e.path)),
            }
        }
        for e in self.entries.iter().filter(|e| e.role == NESTED_OUTPUT && path_inside(&e.path)) {
            if e.file_name() != "summary.json" || e.dir().is_empty() || e.format != "json" || e.scale != "none" {
                out.push(format!("{}: a nested output is listed by its own summary.json, in a subdirectory, as json with scale none", e.path));
                continue;
            }
            if self.entries.iter().any(|x| x.path != e.path && x.path.starts_with(&format!("{}/", e.dir()))) {
                out.push(format!("{}/: a nested output's files are listed by its own manifest, not by this one", e.dir()));
            }
            if depth >= MAX_NESTING {
                out.push(format!("{}: nested outputs more than {MAX_NESTING} levels deep are not checked", e.path));
                continue;
            }
            let sub = dir.join(e.dir());
            match Manifest::read(&sub) {
                Ok(m) => out.extend(m.check_at(&sub, verify_digests, depth + 1)?.into_iter().map(|p| format!("{}/: {p}", e.dir()))),
                Err(x) => out.push(format!("{}: not a scrust output's manifest: {x}", e.path)),
            }
        }
        out.extend(naming_problems(&self.entries));
        Ok(out)
    }

    pub fn to_json(&self) -> Json {
        Json::obj(vec![
            ("schema", Json::Str(SCHEMA.into())),
            ("version", Json::Int(VERSION)),
            ("documentation", Json::Str(DOCUMENTATION.into())),
            ("command", Json::Str(self.command.clone())),
            ("provenance", self.provenance.to_json()),
            ("entries", Json::List(self.entries.iter().map(Entry::to_json).collect())),
        ])
    }

    /// Reads a manifest (`summary.json`'s `files` section).
    pub fn from_json(j: &Json) -> Result<Manifest, String> {
        if j.get("schema").and_then(Json::as_str) != Some(SCHEMA) {
            return Err(format!("not a {SCHEMA} manifest"));
        }
        let v = j.get("version").and_then(Json::as_f64).ok_or("no version")? as i64;
        if v != VERSION {
            return Err(format!("{SCHEMA} version {v}; this build reads version {VERSION}"));
        }
        let p = j.get("provenance").ok_or("no provenance")?;
        let strs = |x: Option<&Json>| -> Vec<String> {
            match x {
                Some(Json::Vector(v) | Json::List(v)) => v.iter().filter_map(|s| s.as_str().map(String::from)).collect(),
                _ => Vec::new(),
            }
        };
        let inputs = match p.get("inputs") {
            Some(Json::List(v) | Json::Vector(v)) => v.iter().map(Input::from_json).collect::<Result<Vec<_>, _>>()?,
            _ => Vec::new(),
        };
        let provenance = Provenance {
            scrust_version: p.get("scrust_version").and_then(Json::as_str).unwrap_or("").into(),
            command_line: strs(p.get("command_line")),
            command: p.get("command").and_then(Json::as_str).unwrap_or("").into(),
            inputs,
            reproduces: strs(p.get("reproduces")),
            parameters: p.get("parameters").cloned().unwrap_or(Json::Null),
            seeds: p.get("seeds").cloned().unwrap_or(Json::Null),
        };
        let entries = match j.get("entries") {
            Some(Json::List(v) | Json::Vector(v)) => v.iter().map(Entry::from_json).collect::<Result<Vec<_>, _>>()?,
            _ => return Err("no entries".into()),
        };
        Ok(Manifest { command: j.get("command").and_then(Json::as_str).unwrap_or("").into(), provenance, entries })
    }

    /// The manifest of the output directory `dir`: its `summary.json`'s `files` section.
    pub fn read(dir: &Path) -> Result<Manifest, String> {
        let p = dir.join("summary.json");
        let text = std::fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()))?;
        let j = Json::parse(&text).map_err(|e| format!("{}: {e}", p.display()))?;
        Manifest::from_json(j.get("files").ok_or_else(|| format!("{}: no files section", p.display()))?)
    }

    pub fn entry(&self, path: &str) -> Option<&Entry> {
        self.entries.iter().find(|e| e.path == path)
    }
}

/// How a file of `format` loads with standard readers, in R and in Python.
fn how_to_load(e: &Entry) -> (String, String) {
    let p = &e.path;
    match (e.format.as_str(), e.header) {
        ("matrix_market", _) => (format!("Matrix::readMM(\"{p}\")"), format!("scipy.io.mmread(\"{p}\")")),
        ("tsv", Some(true)) => (format!("read.delim(\"{p}\", check.names = FALSE)"), format!("pandas.read_csv(\"{p}\", sep=\"\\t\")")),
        ("tsv", _) => {
            if e.shape.as_ref().is_some_and(|s| s.get(1) == Some(&1)) {
                (format!("readLines(\"{p}\")"), format!("pandas.read_csv(\"{p}\", sep=\"\\t\", header=None)[0]"))
            } else {
                (format!("read.delim(\"{p}\", header = FALSE)"), format!("pandas.read_csv(\"{p}\", sep=\"\\t\", header=None)"))
            }
        }
        ("npy", _) => (format!("RcppCNPy::npyLoad(\"{p}\") (or reticulate's numpy)"), format!("numpy.load(\"{p}\")")),
        ("json", _) => (format!("jsonlite::fromJSON(\"{p}\")"), format!("json.load(open(\"{p}\"))")),
        ("annoy", _) => (format!("RcppAnnoy: new(AnnoyEuclidean, dims)$load(\"{p}\")"), format!("annoy.AnnoyIndex(dims, \"euclidean\").load(\"{p}\")")),
        _ => (format!("readLines(\"{p}\")"), format!("open(\"{p}\").read()")),
    }
}

/// `README.txt` for an output directory: `title`, `intro` (paragraphs), then every file the
/// manifest lists, with its role, scale, shape, identifiers and how to load it with standard
/// readers.
pub fn readme(m: &Manifest, title: &str, intro: &[String]) -> String {
    let mut t = format!("{title}\n{}\n\n", "=".repeat(title.chars().count()));
    for para in intro {
        t.push_str(para);
        t.push_str("\n\n");
    }
    let mut dirs: Vec<&str> = m.entries.iter().map(|e| e.dir()).collect();
    dirs.dedup();
    let mut seen: BTreeSet<&str> = BTreeSet::new();
    for d in dirs {
        if !seen.insert(d) {
            continue;
        }
        t.push_str(&format!("{}\n{}\n\n", if d.is_empty() { "./" } else { d }, "-".repeat(if d.is_empty() { 2 } else { d.chars().count() })));
        for e in m.entries.iter().filter(|e| e.dir() == d) {
            t.push_str(&format!("{}\n    {}\n", e.path, e.description));
            let mut facts = vec![format!("role {}", e.role)];
            if e.scale != "none" {
                facts.push(format!("scale {}", e.scale));
            }
            if let Some(s) = &e.shape {
                facts.push(format!("shape {}", s.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" x ")));
            }
            if let Some(n) = e.nnz {
                facts.push(format!("{n} stored entries"));
            }
            let ids = |what: &str, i: &Ids| {
                let col = match &i.column {
                    Column::Name(n) => format!("column {n}"),
                    Column::Number(k) => format!("column {k}"),
                };
                format!("{what}: {} ({}, {col})", i.kind, i.file)
            };
            if let Some(r) = &e.rows {
                facts.push(ids("rows", r));
            }
            if let Some(c) = &e.columns {
                facts.push(ids("columns", c));
            }
            t.push_str(&format!("    {}\n", facts.join("; ")));
            if e.path != "README.txt" && e.role != "stage_dump" {
                let (r, py) = how_to_load(e);
                t.push_str(&format!("    R: {r}\n    Python: {py}\n"));
            }
            t.push('\n');
        }
    }
    t
}

/// The command line as the normalised matrices' comment lines and the manifest give it: `scrust`,
/// the subcommand and its arguments, without the options in `drop` (each with its value, unless
/// it is one of `switches`), which change no value of any file the command writes.
pub fn canonical_command(subcommand: &str, args: &[String], drop: &[&str], switches: &[&str]) -> String {
    let mut out = vec!["scrust".to_string(), subcommand.to_string()];
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        let name = a.strip_prefix("--").unwrap_or("");
        let takes_value = !switches.contains(&name);
        if drop.contains(&name) {
            i += if takes_value { 2 } else { 1 };
            continue;
        }
        out.push(a.clone());
        if takes_value && !name.is_empty() {
            if let Some(v) = args.get(i + 1) {
                out.push(v.clone());
            }
            i += 2;
        } else {
            i += 1;
        }
    }
    out.iter().map(|w| if w.is_empty() || w.contains(char::is_whitespace) { format!("'{w}'") } else { w.clone() }).collect::<Vec<_>>().join(" ")
}

/// Entries for a stage-dump directory (`rel`, e.g. `stages`, under `out`): every file, with the
/// `.npy` arrays' shapes from their headers; R's layout is described by `stages/index.tsv`.
pub fn stage_entries(out: &Path, rel: &str) -> io::Result<Vec<Entry>> {
    let mut v = Vec::new();
    let dir = out.join(rel);
    if !dir.is_dir() {
        return Ok(v);
    }
    for f in list_files(&dir)? {
        let path = format!("{rel}/{f}");
        let p = dir.join(&f);
        let e = if f.ends_with(".npy") {
            let shape = crate::npy::Npy::read_header(&p).map_err(|e| io::Error::other(e.to_string()))?.1;
            let what = if f.contains("_partitions") { "every random start's partition (starts x cells), as scrust cluster --keep-stages keeps them" } else { "a stage dump in the instrumented R runs' layout (stages/index.tsv says what each block holds)" };
            Entry::new(&path, if f.contains("_partitions") { "start_partitions" } else { "stage_dump" }, "none", "npy", what).shape(&shape)
        } else if f == "index.tsv" || f == "rng.tsv" {
            Entry::new(&path, "stage_index", "none", "tsv", "the stage dumps' index (which array, context, rows and stage)").header(true)
        } else if f == "stages.json" {
            Entry::new(&path, "stage_scalars", "none", "json", "the stage dumps' scalars, per context")
        } else if f.ends_with(".ann") {
            Entry::new(&path, "stage_dump", "none", "annoy", "Annoy's index as save() writes it (stage K4)")
        } else {
            Entry::new(&path, "stage_dump", "none", "text", "a stage dump")
        };
        v.push(e);
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn provenance() -> Provenance {
        Provenance {
            scrust_version: "0.1.0".into(),
            command_line: vec!["scrust".into(), "x".into()],
            command: "scrust x".into(),
            inputs: vec![Input { option: "counts".into(), path: "in".into(), files: vec![InputFile { name: "matrix.mtx.gz".into(), bytes: 3, sha256: "ab".into() }] }],
            reproduces: vec!["R".into()],
            parameters: Json::obj(vec![("k", Json::Int(1))]),
            seeds: Json::obj(vec![("seed", Json::Int(2))]),
        }
    }

    #[test]
    fn the_naming_rules() {
        let raw = |p: &str, role: &str| Entry::new(p, role, if role == "raw_counts" { "raw_counts" } else { "none" }, "tsv", "");
        let ok = vec![
            raw("counts/matrix.mtx.gz", "raw_counts"),
            raw("counts/barcodes.tsv.gz", "raw_counts_cells"),
            raw("counts/features.tsv.gz", "raw_counts_features"),
            Entry::new("sct/corrected_counts.mtx.gz", "sct_counts", "corrected_counts", "matrix_market", ""),
            Entry::new("sct/sct_features.tsv.gz", "sct_features", "none", "tsv", ""),
            Entry::new("metadata.tsv", "cell_metadata", "none", "tsv", ""),
        ];
        assert!(naming_problems(&ok).is_empty(), "{:?}", naming_problems(&ok));
        // a normalised matrix under a 10x name, or beside one, or ending in matrix.mtx.gz
        for (bad, want) in [
            (Entry::new("sct/matrix.mtx.gz", "sct_counts", "corrected_counts", "matrix_market", ""), 3),
            (Entry::new("counts/log1p.mtx.gz", "sct_data", "log1p_corrected_counts", "matrix_market", ""), 1),
            (Entry::new("sct/corrected_counts.matrix.mtx.gz", "sct_counts", "corrected_counts", "matrix_market", ""), 1),
            (Entry::new("barcodes.tsv.gz", "sct_cells", "none", "tsv", ""), 2 + 1),
        ] {
            let mut v = ok.clone();
            v.push(bad.clone());
            let p = naming_problems(&v);
            assert!(p.len() >= want.min(1), "{bad:?}: {p:?}");
            assert!(p.iter().all(|m| m.contains(&bad.path) || m.contains("10x")), "{p:?}");
        }
    }

    #[test]
    fn round_trip_and_checks() {
        let dir = crate::test_dir("manifest");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("sct")).unwrap();
        std::fs::write(dir.join("metadata.tsv"), "barcode\tx\nA\t1\n").unwrap();
        std::fs::write(dir.join("sct/sct_cells.tsv.gz"), b"zz").unwrap();
        std::fs::write(dir.join("summary.json"), "{}").unwrap();
        let mut m = Manifest::new("cluster", provenance());
        m.push(Manifest::summary_entry());
        m.push(Entry::new("metadata.tsv", "cell_metadata", "none", "tsv", "per cell").header(true).shape(&[1, 2]).rows(Ids::named("cells", "metadata.tsv", "barcode")));
        m.push(Entry::new("sct/sct_cells.tsv.gz", "sct_cells", "none", "tsv", "cells").header(false).shape(&[1, 1]));
        m.fill_digests(&dir, 3).unwrap();
        assert_eq!(m.entries[1].sha256.as_deref(), Some(crate::sha256::hex(b"barcode\tx\nA\t1\n").as_str()));
        assert!(m.check(&dir, true).unwrap().is_empty(), "{:?}", m.check(&dir, true).unwrap());
        let text = crate::json::to_file_text(&Json::obj(vec![("files", m.to_json())]));
        let back = Manifest::from_json(Json::parse(&text).unwrap().get("files").unwrap()).unwrap();
        assert_eq!(back, m);
        // the same digests at one thread
        let mut m1 = m.clone();
        m1.entries.iter_mut().for_each(|e| e.sha256 = None);
        m1.fill_digests(&dir, 1).unwrap();
        assert_eq!(m1, m);
        // an unlisted file, a missing one, a changed one
        std::fs::write(dir.join("extra.txt"), "x").unwrap();
        std::fs::write(dir.join("metadata.tsv"), "barcode\tx\nA\t2\n").unwrap();
        m.push(Entry::new("gone.tsv", "x", "none", "tsv", ""));
        let p = m.check(&dir, true).unwrap();
        assert!(p.iter().any(|s| s.starts_with("extra.txt")) && p.iter().any(|s| s.starts_with("gone.tsv")) && p.iter().any(|s| s.starts_with("metadata.tsv")), "{p:?}");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// A nested output (M15): the parent lists only its `summary.json`; the check skips its files
    /// and checks its own manifest instead; given digests are kept; no listed path, nested or
    /// not, may leave the output.
    #[test]
    fn nested_outputs() {
        let dir = crate::test_dir("manifest_nested");
        let _ = std::fs::remove_dir_all(&dir);
        let sub = dir.join("singlets_clusters");
        std::fs::create_dir_all(sub.join("counts")).unwrap();
        std::fs::write(sub.join("metadata.tsv"), "barcode\tx\nA\t1\n").unwrap();
        std::fs::write(sub.join("counts/matrix.mtx.gz"), b"m").unwrap();
        let mut inner = Manifest::new("cluster", provenance());
        inner.push(Manifest::summary_entry());
        inner.push(Entry::new("metadata.tsv", "cell_metadata", "none", "tsv", "per cell").header(true).shape(&[1, 2]));
        inner.push(Entry::new("counts/matrix.mtx.gz", "raw_counts", "raw_counts", "matrix_market", "counts"));
        inner.fill_digests(&sub, 1).unwrap();
        let write_summary = |d: &Path, m: &Manifest| std::fs::write(d.join("summary.json"), crate::json::to_file_text(&Json::obj(vec![("files", m.to_json())]))).unwrap();
        write_summary(&sub, &inner);
        std::fs::write(dir.join("cells.tsv"), "barcode\npANN\n").unwrap();
        let mut m = Manifest::new("doubletfinder", provenance());
        m.push(Manifest::summary_entry());
        let mut cells = Entry::new("cells.tsv", "doublet_calls", "none", "tsv", "calls").header(true);
        // a digest given by the writer is kept, not recomputed
        cells.bytes = Some(15);
        cells.sha256 = Some(crate::sha256::hex(b"barcode\npANN\n"));
        m.push(cells);
        m.push(Entry::new("singlets_clusters/summary.json", NESTED_OUTPUT, "none", "json", "nested"));
        m.fill_digests(&dir, 2).unwrap();
        assert_eq!(m.entries[1].bytes, Some(15), "the given size is kept");
        assert!(m.entries[2].sha256.is_some(), "the nested summary.json is digested");
        m.entries[1].bytes = Some(13);
        write_summary(&dir, &m);
        let p = m.check(&dir, true).unwrap();
        assert!(p.is_empty(), "{p:?}");
        // the nested output's own problems are the parent's, under its directory, digests too
        std::fs::write(sub.join("stray.txt"), "x").unwrap();
        let p = m.check(&dir, false).unwrap();
        assert_eq!(p, vec!["singlets_clusters/: stray.txt: in the directory, not in the manifest".to_string()]);
        std::fs::remove_file(sub.join("stray.txt")).unwrap();
        std::fs::write(sub.join("metadata.tsv"), "barcode\tx\nA\t2\n").unwrap();
        let p = m.check(&dir, true).unwrap();
        assert!(p.len() == 1 && p[0].starts_with("singlets_clusters/: metadata.tsv: 14 bytes, sha256 ") && p[0].contains("the manifest says 14 bytes"), "{p:?}");
        std::fs::write(sub.join("metadata.tsv"), "barcode\tx\nA\t1\n").unwrap();
        // a nested output's file listed by the parent too, and a nested entry not named summary.json
        let mut m2 = m.clone();
        m2.push(Entry::new("singlets_clusters/metadata.tsv", "cell_metadata", "none", "tsv", ""));
        assert!(m2.check(&dir, false).unwrap().iter().any(|s| s.contains("listed by its own manifest")));
        let mut m3 = m.clone();
        m3.entries[2].path = "singlets_clusters/metadata.tsv".into();
        assert!(m3.check(&dir, false).unwrap().iter().any(|s| s.contains("by its own summary.json")));
        // a nested output can't escape its parent: `..` and absolute paths are refused, not followed
        let outside = dir.with_file_name("manifest_nested_outside");
        let _ = std::fs::remove_dir_all(&outside);
        std::fs::create_dir_all(&outside).unwrap();
        write_summary(&outside, &Manifest::new("cluster", provenance()));
        for bad in ["../manifest_nested_outside/summary.json".to_string(), format!("{}/summary.json", outside.display()), "singlets_clusters/./summary.json".into(), "singlets_clusters//summary.json".into()] {
            let mut m4 = m.clone();
            m4.entries[2].path = bad.clone();
            let p = m4.check(&dir, false).unwrap();
            assert!(p.iter().any(|s| s == &format!("{bad}: not a relative path inside the output (absolute, or with a `..`, `.` or empty component)")), "{bad}: {p:?}");
            assert!(!p.iter().any(|s| s.contains("not a scrust output's manifest") || s.contains("listed, not in")), "{bad}: not followed: {p:?}");
        }
        // nor may any other file
        let mut m5 = m.clone();
        m5.push(Entry::new("../manifest_nested_outside/summary.json", "cell_metadata", "none", "json", ""));
        assert!(m5.check(&dir, false).unwrap().iter().any(|s| s.starts_with("../manifest_nested_outside/summary.json: not a relative path inside")));
        // a nested entry with another scale
        let mut m6 = m.clone();
        m6.entries[2].scale = "raw_counts".into();
        assert!(m6.check(&dir, false).unwrap().iter().any(|s| s.contains("as json with scale none")));
        assert!(path_inside("a/b.tsv") && path_inside("summary.json") && !path_inside("") && !path_inside("/a") && !path_inside("a/../b") && !path_inside("a/"));
        std::fs::remove_dir_all(&outside).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// Nested outputs are checked level by level, up to [`MAX_NESTING`] levels deep.
    #[test]
    fn nested_outputs_nest() {
        let dir = crate::test_dir("manifest_nest_deep");
        let _ = std::fs::remove_dir_all(&dir);
        // level k holds a.tsv and, below level MAX_NESTING + 1, the nested output n/
        let levels = MAX_NESTING + 2;
        for k in (0..levels).rev() {
            let d = (0..k).fold(dir.clone(), |d, _| d.join("n"));
            std::fs::create_dir_all(&d).unwrap();
            std::fs::write(d.join("a.tsv"), "x\n").unwrap();
            let mut m = Manifest::new("cluster", provenance());
            m.push(Manifest::summary_entry());
            m.push(Entry::new("a.tsv", "cell_metadata", "none", "tsv", "a").header(true));
            if k + 1 < levels {
                m.push(Entry::new("n/summary.json", NESTED_OUTPUT, "none", "json", "nested"));
            }
            m.fill_digests(&d, 1).unwrap();
            std::fs::write(d.join("summary.json"), crate::json::to_file_text(&Json::obj(vec![("files", m.to_json())]))).unwrap();
        }
        let m = Manifest::read(&dir).unwrap();
        let p = m.check(&dir, true).unwrap();
        let deep = format!("{}n/summary.json: nested outputs more than {MAX_NESTING} levels deep are not checked", "n/: ".repeat(MAX_NESTING));
        assert_eq!(p, vec![deep]);
        // a problem at level 3 is reported under its directories
        std::fs::write(dir.join("n/n/n/a.tsv"), "y\n").unwrap();
        let p = m.check(&dir, true).unwrap();
        assert!(p.iter().any(|s| s.starts_with("n/: n/: n/: a.tsv: 2 bytes, sha256 ")), "{p:?}");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn the_canonical_command() {
        let a: Vec<String> = ["--counts", "in dir", "--out", "o", "--threads", "8", "--export", "--preset", "workflow", "--keep-stages"].iter().map(|s| s.to_string()).collect();
        assert_eq!(canonical_command("cluster", &a, &["out", "threads", "keep-stages"], &["export", "keep-stages"]), "scrust cluster --counts 'in dir' --export --preset workflow");
    }
}
