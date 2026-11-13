#!/bin/bash
# M8's R-vs-R floor against the saved objects (validation/R/floor_clustering.R): every portable run
# against its saved object, output by output, into <out.tsv>, and the report's tables into
# <tables.md>. Reads only plain files, but the tier-P scale.data and PCA arrays need a job (16 GB
# cover them). It needs floor_scale_data.sbatch's measurements of the objects whose scale.data the
# reference runs kept only as a digest (dump level light), and no particular CPU type, so
# --partition=shared,sapphire starts it sooner.
# From the repository's top directory:
#   sbatch validation/slurm/floor_clustering.sbatch <out.tsv> <tables.md>
#
#SBATCH --job-name=scrust-cl-floor-table
#SBATCH -A informatics
#SBATCH --partition=shared,sapphire
#SBATCH --exclude=holy8a32608,holy8a24102
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0-04:00:00
#SBATCH --output=%x-%j.sbatch.log

set -euo pipefail

source /n/holylfs05/LABS/informatics/Lab/projects/afreedman/scRNAseq-preprocessing-analysis/results/_regen_mt/node_guard.sh
cd "${SCRUST_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ $# -eq 2 ] || { echo "usage: floor_clustering.sbatch <out.tsv> <tables.md>" >&2; exit 2; }
R_ENV=/n/holylfs05/LABS/informatics/Users/afreedman/smk-conda-envs/d3ca7f730a70e6f748bb65a12a5ed3f9_
env -u R_LIBS_USER -u R_LIBS_SITE -u R_LIBS -u PERL5LIB PYTHONNOUSERSITE=1 \
    /usr/bin/time -v -o "$2.time_v.txt" "$R_ENV/bin/Rscript" validation/R/floor_clustering.R "$1" > "$2.partial"
mv "$2.partial" "$2"
