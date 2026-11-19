#!/usr/bin/env python3
"""Make the benchmark copies of the preprocessing workflow's R scripts (scRust data/benchmarks).

Each copy is the workflow's script, byte for byte, with phase-marker lines inserted before
anchor lines; nothing else changes. The script checks that removing the inserted lines gives
back the source exactly, and records each source's git blob id and sha256.

  PYTHONNOUSERSITE=1 python3 -I make_workflow_copies.py

Sources: <pipeline>/workflow/scripts/<name> at the pipeline repo's HEAD (read only).
Copies:  validation/benchmark/workflow/scripts/<name>, plus SOURCES.tsv.
"""

import hashlib
import os
import subprocess
import sys

PIPE = "/n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq_preprocessing/scRNAseq-preprocessing"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "workflow", "scripts")

HEAD = 'source("workflow/scripts/bench_markers.R")\nbench_phase("packages", "startup")\n'

def chain_markers(var, kind="compute"):
    """The ported chain's markers, for the Seurat object named var."""
    return [
        ('%s[["percent.mt"]] <- PercentageFeatureSet(' % var, 'bench_phase("percent_mt", "%s")' % kind),
        ("%s <- SCTransform(" % var, 'bench_phase("sctransform", "%s")' % kind),
        ("require_min_cells_for_pca(%s" % var, 'bench_phase("pca", "%s")' % kind),
        ("%s <- RunUMAP(" % var, 'bench_phase("umap", "unported")'),
        ("%s <- FindNeighbors(" % var, 'bench_phase("neighbors", "%s")' % kind),
        ("%s <- FindClusters(" % var, 'bench_phase("clusters", "%s")' % kind),
        ("%s <- add_silhouette_to_metadata(" % var, 'bench_phase("silhouette", "%s")' % kind),
    ]


def write_markers(saved_call, meta_call):
    return [
        (saved_call, 'bench_phase("write_rds", "write")'),
        (meta_call, 'bench_phase("write_cluster_ids", "write")'),
    ]


SCRIPTS = {
    "tenx2seuratrds.R": [
        ("filtered_loaded <- Seurat::Read10X(filtered)", 'bench_phase("read", "read")'),
    ] + chain_markers("seurat_obj") + write_markers(
        "saveRDS(seurat_obj,file=output)", "write_cluster_metadata(seurat_obj,"),
    "soupx.R": [
        ("filtered_matrix <- Seurat::Read10X(filtered)", 'bench_phase("read", "read")'),
        # SoupX's other inputs, read only for SoupX: their own read phases
        ("raw_matrix <- Seurat::Read10X(raw)", 'bench_phase("read_raw", "read")'),
        ("seurat_base <- readRDS(seurat_base)", 'bench_phase("read_base_rds", "read")'),
        ("soup_channel <- SoupX::SoupChannel(", 'bench_phase("soupx", "unported")'),
    ] + chain_markers("seurat_soupx") + write_markers(
        "saveRDS(seurat_soupx,file=output)", "write_cluster_metadata(seurat_soupx,"),
    "emptydrops.R": [
        ("droputil_rawdata <-read10xCounts(raw)", 'bench_phase("read", "read")'),
        ("set.seed(WORKFLOW_SEED)", 'bench_phase("emptydrops", "unported")'),
    ] + chain_markers("seurat_droputil_filtered") + write_markers(
        "saveRDS(seurat_droputil_filtered,file=seurat_output)",
        "write_cluster_metadata(seurat_droputil_filtered,") + [
        ("dir.create(matrix_outdir", 'bench_phase("write_matrix", "write")'),
    ],
    "cellbender2seurat.R": [
        ("mat <- Read_CellBender_h5_Mat(cellbender_h5)", 'bench_phase("read", "read")'),
    ] + chain_markers("seurat") + write_markers(
        "saveRDS(seurat,file=output)", "write_cluster_metadata(seurat,"),
    "doubletfinder.R": [
        ("seurat <- readRDS(seurat)", 'bench_phase("read", "read")'),
        ("sweep.res.list <- paramSweep(", 'bench_phase("df_sweep", "compute")'),
        ("sweep.stats <- summarizeSweep(", 'bench_phase("df_pk", "compute")'),
        ("nExp_poi <- round(", 'bench_phase("df_final", "compute")'),
        ("classifier_colname <- paste(", 'bench_phase("df_subset", "compute")'),
    ] + [(a, l.replace('"compute"', '"recluster"')) for a, l in chain_markers("seurat_singlets")] + write_markers(
        "saveRDS(seurat_singlets,file=output)", "write_cluster_metadata(seurat_singlets,"),
    "silhouette_utils.R": [
        ("  purity <- bluster::neighborPurity(", '  bench_subphase("purity")'),
    ],
}

# anchors that occur more than once: which occurrence (1-based) gets the marker
OCCURRENCE = {
    ("emptydrops.R", "set.seed(WORKFLOW_SEED)"): 1,
}


def git_blob(path):
    rel = os.path.relpath(path, PIPE)
    out = subprocess.check_output(
        ["timeout", "60", "git", "--no-optional-locks", "-C", PIPE, "rev-parse", "HEAD:" + rel])
    return out.decode().strip()


def main():
    os.makedirs(OUT, exist_ok=True)
    head = subprocess.check_output(
        ["timeout", "60", "git", "--no-optional-locks", "-C", PIPE, "rev-parse", "HEAD"]).decode().strip()
    dirty = subprocess.check_output(
        ["timeout", "60", "git", "--no-optional-locks", "-C", PIPE, "status", "--porcelain", "--",
         "workflow/scripts", "workflow/rules"]).decode().strip()
    rows = []
    for name, markers in SCRIPTS.items():
        src = os.path.join(PIPE, "workflow", "scripts", name)
        raw = open(src, "rb").read()
        text = raw.decode()
        lines = text.split("\n")
        out = []
        inserted = 0
        used = set()
        seen_count = {}
        for line in lines:
            for anchor, marker in markers:
                if line.startswith(anchor):
                    seen_count[anchor] = seen_count.get(anchor, 0) + 1
                    want = OCCURRENCE.get((name, anchor), 1)
                    if want is None or seen_count[anchor] != want:
                        continue
                    if anchor in used:
                        sys.exit("%s: anchor used twice: %s" % (name, anchor))
                    used.add(anchor)
                    out.append(marker)
                    inserted += 1
            out.append(line)
        missing = [a for a, _ in markers if a not in used]
        if missing:
            sys.exit("%s: anchors not found: %s" % (name, missing))
        body = "\n".join(out)
        if name != "silhouette_utils.R":
            body = HEAD + body
            if not body.endswith("\n"):
                body += "\n"
            body += "bench_end()\n"
        # check: dropping the inserted lines gives the source back, byte for byte
        marker_lines = set(m for _, m in markers) | set(HEAD.splitlines()) | {"bench_end()"}
        kept = [l for l in body.split("\n") if l not in marker_lines]
        back = "\n".join(kept)
        if back.rstrip("\n") != text.rstrip("\n"):
            sys.exit("%s: the copy minus its markers differs from the source" % name)
        dst = os.path.join(OUT, name)
        tmp = dst + ".tmp"
        with open(tmp, "w") as fh:
            fh.write(body)
        os.rename(tmp, dst)  # atomic: an R process that has the old file open keeps reading it
        rows.append((name, src, git_blob(src), hashlib.sha256(raw).hexdigest(),
                     hashlib.sha256(body.encode()).hexdigest(), str(inserted + (0 if name == "silhouette_utils.R" else 3))))
        print("%-22s %2d marker lines" % (name, inserted + (0 if name == "silhouette_utils.R" else 3)))
    with open(os.path.join(HERE, "workflow", "SOURCES.tsv"), "w") as fh:
        fh.write("# pipeline repo %s HEAD %s; workflow/scripts and workflow/rules %s\n" % (
            PIPE, head, "clean" if not dirty else "MODIFIED: " + dirty.replace("\n", "; ")))
        fh.write("copy\tsource\tgit_blob\tsource_sha256\tcopy_sha256\tlines_added\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
