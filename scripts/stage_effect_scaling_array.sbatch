#!/bin/bash
#SBATCH --job-name=motco-ses
#SBATCH --output=logs/motco-ses-%A_%a.out
#SBATCH --error=logs/motco-ses-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Stage-effect sample-size scaling — SLURM array template for the two profiles.
#
# Usage (from a clean worktree at the run's pushed commit, submitted from the worktree's root; every path
# absolute; run `calibrate` first, pre-create each profile's settings.json so no two tasks race to write it, and
# mkdir logs before sbatch — SLURM opens --output before the script runs):
#   mkdir -p logs
#   # the parent's committed stage_span aggregates are read from the worktree's results/; the band settings bind
#   # /abs/run/calibration/calibration.json by sha256: copy the calibration directory the smoke was built from
#   # into the run root (or run `calibrate --out-dir /abs/run --workers 8` there, once, before any band unit)
#   mkdir -p /abs/run && cp -r /local/run/calibration /abs/run/
#   for P in saturated band; do
#     PYTHONPATH=$PWD/src python scripts/stage_effect_scaling_check.py run --profile $P --out-dir /abs/run \
#       --shard 0/1 --limit 0
#   done
#   sbatch --array=0-199%64 --export=ALL,SES_PROFILE=saturated,SES_OUT=/abs/run,N_SHARDS=200 \
#     scripts/stage_effect_scaling_array.sbatch
#   sbatch --array=0-139%64 --export=ALL,SES_PROFILE=band,SES_OUT=/abs/run,N_SHARDS=140 \
#     scripts/stage_effect_scaling_array.sbatch
#   # after every task of both arrays completes:
#   python scripts/stage_effect_scaling_check.py merge --profile saturated --out-dir /abs/run
#   python scripts/stage_effect_scaling_check.py merge --profile band --out-dir /abs/run
#   python scripts/stage_effect_scaling_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   SES_PROFILE   saturated | band.
#   SES_OUT       Absolute run root; shards land in $SES_OUT/<profile>/records/, the band reads
#                 $SES_OUT/calibration/calibration.json.
#   N_SHARDS      Total shards = --array width.
# Optional:
#   SES_ARGS      Extra run flags ("--limit 1" for timing). They must not change the settings between
#                 submissions: a resume into a directory with different settings is refused.
#   SES_DRY_RUN   When set, print the run command instead of executing it.
#
# PYTHONPATH is set here to the submitting worktree's src/, so the run measures that revision's package.
# Each task runs its shard's units one at a time (one unit per CPU, BLAS pinned to one thread): a unit
# generates its dataset, fits production once, evaluates it, then runs the split test, the interaction scan
# (and the false-discovery reading from its F values) and the matched-PCA evaluation at the unit's RRPP seed.
# Shards are round-robin over a design-major order, so each mixes the n = 1,280 units with the cheaper ones, and
# resumable: resubmitting failed array ids skips complete units. The `%64` throttle caps concurrent tasks; change
# it mid-run with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${SES_PROFILE:?SES_PROFILE environment variable is required (saturated or band).}"
: "${SES_OUT:?SES_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

case "$SES_PROFILE" in
    saturated | band) ;;
    *)
        echo "SES_PROFILE must be saturated or band, got '$SES_PROFILE'." >&2
        exit 2
        ;;
esac
case "$SES_OUT" in
    /*) ;;
    *)
        echo "SES_OUT must be an absolute path, got '$SES_OUT'." >&2
        exit 2
        ;;
esac

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$REPO_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$REPO_DIR/src"

# shellcheck disable=SC2206
COMMAND=(python "$REPO_DIR/scripts/stage_effect_scaling_check.py" run --profile "$SES_PROFILE" --out-dir "$SES_OUT"
         --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${SES_ARGS:-})

if [ -n "${SES_DRY_RUN:-}" ]; then
    printf 'PYTHONPATH=%q ' "$PYTHONPATH"
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

mkdir -p "$SES_OUT" logs
if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
"${COMMAND[@]}"
