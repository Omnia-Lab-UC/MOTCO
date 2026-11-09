#!/bin/bash
#SBATCH --job-name=motco-occ
#SBATCH --output=logs/motco-occ-%A_%a.out
#SBATCH --error=logs/motco-occ-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Orthonormal span coordinates against the PLS scores — SLURM array template for both profiles.
#
# Usage (from a clean worktree at the run's commit; every path absolute; pre-create each profile's
# settings.json first so no two tasks race to write it, and mkdir logs before sbatch — SLURM opens --output
# before the script runs):
#   mkdir -p logs
#   # the parents (paper-grade profiles and operating tables, span-split settings and production table,
#   # the comparison's calibration) are read from the worktree's examples/ and results/
#   PYTHONPATH=$PWD/src python scripts/orthonormal_coordinates_check.py run --profile paper_grade \
#     --out-dir /abs/run --shard 0/650 --limit 0
#   sbatch --array=0-649%64 --export=ALL,PYTHONPATH=$PWD/src,ORTHO_PROFILE=paper_grade,ORTHO_OUT=/abs/run,\
#   N_SHARDS=650 scripts/orthonormal_coordinates_array.sbatch
#   PYTHONPATH=$PWD/src python scripts/orthonormal_coordinates_check.py run --profile weak \
#     --out-dir /abs/run --shard 0/20 --limit 0
#   sbatch --array=0-19%20 --export=ALL,PYTHONPATH=$PWD/src,ORTHO_PROFILE=weak,ORTHO_OUT=/abs/run,\
#   N_SHARDS=20 scripts/orthonormal_coordinates_array.sbatch
#   # after every task of both arrays completes:
#   python scripts/orthonormal_coordinates_check.py merge --profile paper_grade --out-dir /abs/run
#   python scripts/orthonormal_coordinates_check.py merge --profile weak --out-dir /abs/run
#   python scripts/orthonormal_coordinates_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   ORTHO_PROFILE  paper_grade | weak.
#   ORTHO_OUT      Run root; shards land in $ORTHO_OUT/<profile>/records/.
#   N_SHARDS       Total shards = --array width.
# Optional:
#   ORTHO_ARGS     Extra run flags ("--limit 1" for timing). They must not change the settings between
#                  submissions: a resume into a directory with different settings is refused.
#
# Each task runs its shard's units one at a time (one unit per CPU, BLAS pinned to one thread): a unit
# generates the parent's dataset, fits the production stage-PLS space once and evaluates it in score and in
# orthonormal span coordinates at the unit's RRPP seed. Shards are round-robin and resumable: resubmitting
# failed array ids skips complete units. The `%64` throttle caps concurrent tasks; raise it mid-run with
# `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${ORTHO_PROFILE:?ORTHO_PROFILE environment variable is required (paper_grade or weak).}"
: "${ORTHO_OUT:?ORTHO_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

case "$ORTHO_PROFILE" in
    paper_grade | weak) ;;
    *)
        echo "ORTHO_PROFILE must be paper_grade or weak, got '$ORTHO_PROFILE'." >&2
        exit 2
        ;;
esac
case "$ORTHO_OUT" in
    /*) ;;
    *)
        echo "ORTHO_OUT must be an absolute path, got '$ORTHO_OUT'." >&2
        exit 2
        ;;
esac

mkdir -p "$ORTHO_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/orthonormal_coordinates_check.py run --profile "$ORTHO_PROFILE" \
    --out-dir "$ORTHO_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${ORTHO_ARGS:-}
