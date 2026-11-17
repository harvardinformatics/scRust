#!/bin/bash
# Milestone M12: the tables of docs/validation/clustering.md from the validation runs' plain files
# (no R, no Rust): every `scrust validate cluster` report of a label (m12_cluster.sbatch) with the
# thread checks and the purity certificate, or the benchmark (m12_benchmark.sbatch). M11's
# m11_report.sh with K9 and the certificate. From the repository's top directory:
#   bash validation/slurm/m12_report.sh <label> > docs/validation/clustering_m12.tsv
#   bash validation/slurm/m12_report.sh <label> bench > docs/validation/clustering_m12_bench.tsv
# The comparisons: one row per comparison: object, tier, profile, reference run, kind (same: the
# run's own profile; other: the deliberately mismatched kernel), pass, mode, and the checks per
# status (identical, within floor, differing, missing, n/a, flagged and reported), overall and
# per stage group (K1 SCTransform, K2 RunPCA, K3-K6 neighbours, K7 FindClusters, K8 silhouette and
# files, K9 purity, saved, teacher), the run's cells, clusters, wall time and peak RSS, its thread
# check, and its certificate (flagged pairs: kind, step, cells, ulps below the radius squared).
# The benchmark: one row per run, for each object with a bench/ directory: M8's direct production
# R runs without and with RunUMAP (runs/direct-production-sapphirerapids-noumap/ and
# runs/direct-production-sapphirerapids/), then Rust's at 1 and at many threads (bench/t<n>/):
# the whole process's wall time and peak RSS (time_v.txt), and each stage's time from summary.json
# (run.phases: the time since the previous phase; R's chain_end and Rust's missing stages empty).
set -euo pipefail
label=${1:?label}
if [ "${2:-}" = bench ]; then
    # run.phases of a summary.json, as name=seconds since the previous phase, in order
    phases() {
        awk '/^    "phases": \{/{p=1; next} p && /^    \}/{exit}
             p && /^      "[^"]+": \{/{n=$1; gsub(/[":{]/, "", n)}
             p && /"elapsed_s": /{v=$2; gsub(/,/, "", v); printf "%s=%.3f\n", n, v - prev; prev = v}' "$1"
    }
    get() { printf '%s\n' "$1" | sed -n "s/^$2=//p"; }
    wall() { awk -F': ' '/Elapsed \(wall clock\)/{n = split($2, a, ":"); s = 0; for (i = 1; i <= n; i++) s = s * 60 + a[i]; printf "%.2f", s}' "$1"; }
    rss() { awk -F': ' '/Maximum resident set size/{print $2}' "$1"; }
    printf 'object\tside\tthreads\tumap\twall_s\tpeak_rss_kib\tload_s\tsctransform_s\tpca_s\tumap_s\tneighbors_s\tcluster_s\tsilhouette_s\tpurity_s\twrite_s\n'
    for b in data/clustering/*/rust/"$label"/bench; do
        [ -d "$b" ] || continue
        o=$(echo "$b" | cut -d/ -f3)
        for u in no yes; do
            r=data/clustering/$o/runs/direct-production-sapphirerapids
            [ "$u" = no ] && r=$r-noumap
            p=$(phases "$r/summary.json")
            printf '%s\tR\t1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$o" "$u" "$(wall "$r/time_v.txt")" "$(rss "$r/time_v.txt")" \
                "$(get "$p" load)" "$(get "$p" SCTransform)" "$(get "$p" RunPCA)" "$(get "$p" RunUMAP)" "$(get "$p" FindNeighbors)" \
                "$(get "$p" FindClusters)" "$(get "$p" silhouette)" "$(get "$p" purity)" "$(get "$p" write)"
        done
        for t in $(cd "$b" && ls -d t[0-9]* | tr -d t | sort -n); do
            p=$(phases "$b/t$t/summary.json")
            printf '%s\tRust\t%s\t-\t%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\n' "$o" "$t" "$(wall "$b/time_v_t$t.txt")" "$(rss "$b/time_v_t$t.txt")" \
                "$(get "$p" load)" "$(get "$p" sctransform)" "$(get "$p" pca)" "$(get "$p" neighbors)" "$(get "$p" cluster)" \
                "$(get "$p" silhouette)" "$(get "$p" purity)" "$(get "$p" write)"
        done
    done
    exit 0
fi
jget() { grep -m1 "\"$2\": " "$1" | sed -E 's/^[^:]*: "?([^",]*)"?,?$/\1/'; }
# statuses per stage of a report: identical within_floor differ missing na reported, for stages
# matching $2 (the checks of the "checks" list; by_stage's entries carry no status)
count() {
    awk -v pat="$2" '
        /"stage": /{st=$0; sub(/.*"stage": "/, "", st); sub(/".*/, "", st)}
        /"status": /{s=$0; sub(/.*"status": "/, "", s); sub(/".*/, "", s); if (st ~ pat) c[s]++}
        END {printf "%d\t%d\t%d\t%d\t%d\t%d", c["identical"], c["within floor"], c["differ"], c["missing"], c["n/a"], c["flagged, reported"]}' "$1"
}
# the run's certificate: flagged pairs as kind:step:cell-other@ulps, or "-"
certificate() {
    python3 -I -c 'import json, sys
p = json.load(open(sys.argv[1]))["result"].get("purity") or {}
c = p.get("certificate") or {}
print(";".join("%s:%s:%d-%d@%.0f" % (q["kind"], q["step"], q["cell"], q["other"], -q["ulps_from_threshold"]) for q in c.get("pairs", [])) or "-")' "$1"
}
printf 'object\ttier\tprofile\treference\tkind\tpass\tmode'
for g in all K1 K2 K3-K6 K7 K8 K9 saved teacher; do printf '\t%s_identical\t%s_within\t%s_differ\t%s_missing\t%s_na\t%s_reported' "$g" "$g" "$g" "$g" "$g" "$g"; done
printf '\tcells\tclusters\twall_s\tpeak_rss_kib\tthreads\tcertificate\n'
for f in data/clustering/*/rust/"$label"/{portable,production}/vs_*.json; do
    [ -f "$f" ] || continue
    dir=$(dirname "$f")
    profile=$(basename "$dir")
    object=$(echo "$f" | cut -d/ -f3)
    tier=$(awk -F'\t' -v o="$object" '$1 == o {print $2}' validation/clustering_datasets.tsv)
    name=$(basename "$f" .json)
    kind=same
    case $name in vs_other_*) kind=other; ref=${name#vs_other_} ;; *) ref=${name#vs_} ;; esac
    s=$dir/run/summary.json
    cells=$(jget "$s" n_cells)
    clusters=$(jget "$s" n_clusters)
    wall=$(awk '/"run": /{r=1} r && /"elapsed_s": /{gsub(/[",]/,""); print $2; exit}' "$s")
    rss=$(awk '/"run": /{r=1} r && /^    "peak_rss_kib": /{gsub(/[",]/,""); print $2; exit}' "$s")
    threads=$(cat "$dir/threads.txt" 2>/dev/null | sed -E 's/^threads ([0-9]+) and 1: ([0-9]+) files, ([0-9]+) differing.*/\1 and 1: \3 of \2 differ/' || true)
    mode=$(awk '/"grading": \{/{g=1} g && /"mode": /{gsub(/[",]/,""); print $2; exit}' "$f")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s' "$object" "$tier" "$profile" "$ref" "$kind" "$(jget "$f" pass)" "$mode"
    for pat in '.' '^K1$' '^K2$' '^K[3456]$' '^(K7|log)$' '^(K8|files)$' '^K9$' '^saved$' '^teacher$'; do printf '\t%s' "$(count "$f" "$pat")"; done
    printf '\t%s\t%s\t%s\t%s\t%s\t%s\n' "$cells" "$clusters" "$wall" "$rss" "${threads:--}" "$(certificate "$s")"
done
