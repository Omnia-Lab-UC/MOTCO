#!/bin/bash
#SBATCH --job-name=motco-ud
#SBATCH --output=logs/motco-ud-%A_%a.out
#SBATCH --error=logs/motco-ud-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Unbalance at equal dispersion — SLURM array template (one run, no profile).
#
# Usage (from a clean worktree at the run's pushed commit, submitted from the worktree's root; every path
# absolute; pre-create settings.json so no two tasks race to write it, and mkdir logs before sbatch — SLURM
# opens --output before the script runs):
#   mkdir -p logs
#   # the group-dispersion run's committed columns tables and decision are read from the worktree's results/
#   # and bound by sha256
#   PYTHONPATH=$PWD/src python scripts/unbalance_diagnosis_check.py run --out-dir /abs/run --shard 0/1 --limit 0
#   sbatch --array=0-149%48 --export=ALL,UD_OUT=/abs/run,N_SHARDS=150 scripts/unbalance_diagnosis_array.sbatch
#   # after every task completes:
#   python scripts/unbalance_diagnosis_check.py merge --out-dir /abs/run
#   python scripts/unbalance_diagnosis_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   UD_OUT        Absolute run root; shards land in $UD_OUT/records/.
#   N_SHARDS      Total shards = --array width.
# Optional:
#   UD_ARGS       Extra run flags ("--limit 1" for timing). They must not change the settings between
#                 submissions: a resume into a directory with different settings is refused.
#   UD_DRY_RUN    When set, print the run command instead of executing it.
#
# PYTHONPATH is set here to the submitting worktree's src/, so the run measures that revision's package.
# Each task runs its shard's units one at a time (one unit per CPU, BLAS pinned to one thread): a unit
# generates its dataset, fits production once, evaluates it, then runs the split test, the interaction scan
# (and the false-discovery reading from its F values), the matched-PCA evaluation, the within-group draws on a
# refit at the recorded rank and last the df and leverage corrected draws on another refit at that rank, all at
# the unit's RRPP seed. Shards are round-robin over a design-major order, so each mixes the n = 200 units with
# the cheaper ones, and resumable: resubmitting failed array ids skips complete units. The `%48` throttle caps
# concurrent tasks; change it mid-run with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${UD_OUT:?UD_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

case "$UD_OUT" in
    /*) ;;
    *)
        echo "UD_OUT must be an absolute path, got '$UD_OUT'." >&2
        exit 2
        ;;
esac

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$REPO_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$REPO_DIR/src"

# shellcheck disable=SC2206
COMMAND=(python "$REPO_DIR/scripts/unbalance_diagnosis_check.py" run --out-dir "$UD_OUT"
         --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${UD_ARGS:-})

if [ -n "${UD_DRY_RUN:-}" ]; then
    printf 'PYTHONPATH=%q ' "$PYTHONPATH"
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

mkdir -p "$UD_OUT" logs
if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
"${COMMAND[@]}"
