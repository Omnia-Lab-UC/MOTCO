#!/bin/bash
#SBATCH --job-name=motco-pse
#SBATCH --output=logs/motco-pse-%A_%a.out
#SBATCH --error=logs/motco-pse-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Predictive-subspace size-excess diagnosis (fourteen arms per plasmode unit: the four reference arms, the
# within-stage and out-of-fold checks, and the descriptive arms) — SLURM array template.
#
# Usage (from a clean worktree at the run's pushed commit; every path absolute; the export copied to the cluster,
# its file hashes equal to the plasmode's settings'; mkdir logs before sbatch — SLURM opens --output before the
# script runs). The pre-flight writes settings.json and preflight.json once, so no two tasks race to write them,
# and every task refuses to start without a passed pre-flight:
#   mkdir -p logs
#   PYTHONPATH=$PWD/src python scripts/predictive_size_excess_check.py preflight \
#     --export-dir /abs/cohort_30 --out-dir /abs/run
#   sbatch --array=0-99%50 --export=ALL,PSE_OUT=/abs/run,PSE_EXPORT=/abs/cohort_30,N_SHARDS=100 \
#     scripts/predictive_size_excess_array.sbatch
#   # after every task completes:
#   python scripts/predictive_size_excess_check.py merge --out-dir /abs/run
#   python scripts/predictive_size_excess_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   PSE_OUT       Run root (absolute); shards land in $PSE_OUT/records/.
#   PSE_EXPORT    The sead_sexdiff cohort_30 export directory (absolute).
#   N_SHARDS      Total shards = --array width.
# Optional:
#   PSE_REPO      The worktree to run from (absolute; default the submission directory). Its src/ goes first on
#                 PYTHONPATH, so the run measures that worktree's package, not an installed one.
#   PSE_ARGS      Extra run flags ("--limit 1" for timing). They must not change the settings between submissions:
#                 a resume into a directory with different settings is refused.
#   PSE_DRY_RUN   When set, print the run command (and the PYTHONPATH it would use) instead of executing it.
#
# Each task runs its shard's units one at a time (one CPU, BLAS pinned to one thread): a unit builds its labels
# and plant once, fits production once (an anchor unit at the fidelity rank, a plant unit through one CV) and
# evaluates every arm serially at the unit's RRPP seed. Shards are round-robin and resumable: resubmitting failed
# array ids skips complete units. The `%50` throttle caps concurrent tasks; raise it mid-run with
# `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${PSE_OUT:?PSE_OUT environment variable is required.}"
: "${PSE_EXPORT:?PSE_EXPORT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

REPO="${PSE_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"

for path_var in PSE_OUT PSE_EXPORT REPO; do
    case "${!path_var}" in
        /*) ;;
        *)
            echo "$path_var must be an absolute path, got '${!path_var}'." >&2
            exit 2
            ;;
    esac
done

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

# shellcheck disable=SC2206
COMMAND=(python "$REPO/scripts/predictive_size_excess_check.py" run --export-dir "$PSE_EXPORT"
         --out-dir "$PSE_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${PSE_ARGS:-})

if [ -n "${PSE_DRY_RUN:-}" ]; then
    printf 'PYTHONPATH=%q ' "$PYTHONPATH"
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

cd "$REPO"
mkdir -p "$PSE_OUT" logs
if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
"${COMMAND[@]}"
