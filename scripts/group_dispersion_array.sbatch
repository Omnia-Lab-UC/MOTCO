#!/bin/bash
#SBATCH --job-name=motco-gd
#SBATCH --output=logs/motco-gd-%A_%a.out
#SBATCH --error=logs/motco-gd-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Group-specific dispersion — SLURM array template for the two profiles (core, columns).
#
# Usage (from a clean worktree at the run's pushed commit, submitted from the worktree's root; every path
# absolute; pre-create each profile's settings.json so no two tasks race to write it, and mkdir logs before
# sbatch — SLURM opens --output before the script runs):
#   mkdir -p logs
#   # both parents' committed aggregates (the univariate baseline's stage_span profile and the stage-effect
#   # scaling run's band profile) are read from the worktree's results/ and bound by sha256
#   for P in core columns; do
#     PYTHONPATH=$PWD/src python scripts/group_dispersion_check.py run --profile $P --out-dir /abs/run \
#       --shard 0/1 --limit 0
#   done
#   sbatch --array=0-199%64 --export=ALL,GD_PROFILE=core,GD_OUT=/abs/run,N_SHARDS=200 \
#     scripts/group_dispersion_array.sbatch
#   sbatch --array=0-139%64 --export=ALL,GD_PROFILE=columns,GD_OUT=/abs/run,N_SHARDS=140 \
#     scripts/group_dispersion_array.sbatch
#   # after every task of both arrays completes:
#   python scripts/group_dispersion_check.py merge --profile core --out-dir /abs/run
#   python scripts/group_dispersion_check.py merge --profile columns --out-dir /abs/run
#   python scripts/group_dispersion_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   GD_PROFILE    core | columns.
#   GD_OUT        Absolute run root; shards land in $GD_OUT/<profile>/records/.
#   N_SHARDS      Total shards = --array width.
# Optional:
#   GD_ARGS       Extra run flags ("--limit 1" for timing). They must not change the settings between
#                 submissions: a resume into a directory with different settings is refused.
#   GD_DRY_RUN    When set, print the run command instead of executing it.
#
# PYTHONPATH is set here to the submitting worktree's src/, so the run measures that revision's package.
# Each task runs its shard's units one at a time (one unit per CPU, BLAS pinned to one thread): a unit
# generates its dataset, fits production once, evaluates it, then runs the split test, the interaction scan
# (and the false-discovery reading from its F values) and the matched-PCA evaluation at the unit's RRPP seed,
# and last the within-group draws on a refit at the recorded rank. Shards are round-robin over a design-major
# order, so each mixes the n = 500 and n = 320 units with the cheaper ones, and resumable: resubmitting failed
# array ids skips complete units. The `%64` throttle caps concurrent tasks; change it mid-run with
# `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${GD_PROFILE:?GD_PROFILE environment variable is required (core or columns).}"
: "${GD_OUT:?GD_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

case "$GD_PROFILE" in
    core | columns) ;;
    *)
        echo "GD_PROFILE must be core or columns, got '$GD_PROFILE'." >&2
        exit 2
        ;;
esac
case "$GD_OUT" in
    /*) ;;
    *)
        echo "GD_OUT must be an absolute path, got '$GD_OUT'." >&2
        exit 2
        ;;
esac

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$REPO_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$REPO_DIR/src"

# shellcheck disable=SC2206
COMMAND=(python "$REPO_DIR/scripts/group_dispersion_check.py" run --profile "$GD_PROFILE" --out-dir "$GD_OUT"
         --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${GD_ARGS:-})

if [ -n "${GD_DRY_RUN:-}" ]; then
    printf 'PYTHONPATH=%q ' "$PYTHONPATH"
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

mkdir -p "$GD_OUT" logs
if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
"${COMMAND[@]}"
