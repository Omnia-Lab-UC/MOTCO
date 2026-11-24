#!/bin/bash
#SBATCH --job-name=motco-ubc
#SBATCH --output=logs/motco-ubc-%A_%a.out
#SBATCH --error=logs/motco-ubc-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Univariate interaction baseline — SLURM array template for the three profiles.
#
# Usage (from a clean worktree at the run's pushed commit; every path absolute; pre-create each profile's
# settings.json first so no two tasks race to write it, and mkdir logs before sbatch — SLURM opens --output
# before the script runs):
#   mkdir -p logs
#   # the parents' committed files (the scan validation's simulation settings and comparator table, the
#   # design-guidance slice's two configs and operating tables) are read from the worktree's results/ and examples/
#   for P in parent stage_span surgery; do
#     PYTHONPATH=$PWD/src python scripts/univariate_baseline_check.py run --profile $P --out-dir /abs/run \
#       --shard 0/1 --limit 0
#   done
#   sbatch --array=0-39%64 --export=ALL,PYTHONPATH=$PWD/src,UBC_PROFILE=parent,UBC_OUT=/abs/run,N_SHARDS=40 \
#     scripts/univariate_baseline_array.sbatch
#   sbatch --array=0-99%64 --export=ALL,PYTHONPATH=$PWD/src,UBC_PROFILE=stage_span,UBC_OUT=/abs/run,N_SHARDS=100 \
#     scripts/univariate_baseline_array.sbatch
#   sbatch --array=0-119%64 --export=ALL,PYTHONPATH=$PWD/src,UBC_PROFILE=surgery,UBC_OUT=/abs/run,N_SHARDS=120 \
#     scripts/univariate_baseline_array.sbatch
#   # after every task of the three arrays completes:
#   python scripts/univariate_baseline_check.py merge --profile parent --out-dir /abs/run
#   python scripts/univariate_baseline_check.py merge --profile stage_span --out-dir /abs/run
#   python scripts/univariate_baseline_check.py merge --profile surgery --out-dir /abs/run
#   python scripts/univariate_baseline_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   UBC_PROFILE   parent | stage_span | surgery.
#   UBC_OUT       Run root; shards land in $UBC_OUT/<profile>/records/.
#   N_SHARDS      Total shards = --array width.
# Optional:
#   UBC_ARGS      Extra run flags ("--limit 1" for timing). They must not change the settings between
#                 submissions: a resume into a directory with different settings is refused.
#   UBC_DRY_RUN   When set, print the run command instead of executing it.
#
# Each task runs its shard's units one at a time (one unit per CPU, BLAS pinned to one thread): a unit
# generates its dataset, fits production once, evaluates it, then runs the split test, the interaction scan
# (and the false-discovery reading from its F values) and the matched-PCA evaluation at the unit's RRPP seed
# (`concat` too only when UBC_ARGS names it in --recorded-spaces; it is off by default, design D3). Shards are round-robin and resumable: resubmitting failed array ids skips complete units. The
# `%64` throttle caps concurrent tasks; raise it mid-run with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${UBC_PROFILE:?UBC_PROFILE environment variable is required (parent, stage_span or surgery).}"
: "${UBC_OUT:?UBC_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

case "$UBC_PROFILE" in
    parent | stage_span | surgery) ;;
    *)
        echo "UBC_PROFILE must be parent, stage_span or surgery, got '$UBC_PROFILE'." >&2
        exit 2
        ;;
esac

mkdir -p "$UBC_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2206
COMMAND=(python scripts/univariate_baseline_check.py run --profile "$UBC_PROFILE" --out-dir "$UBC_OUT"
         --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${UBC_ARGS:-})

if [ -n "${UBC_DRY_RUN:-}" ]; then
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi
"${COMMAND[@]}"
