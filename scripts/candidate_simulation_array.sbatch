#!/bin/bash
#SBATCH --job-name=motco-csc
#SBATCH --output=logs/motco-csc-%A_%a.out
#SBATCH --error=logs/motco-csc-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Simulation-only candidate comparison — SLURM array template for the five profiles (committed, conditions,
# dispersion, surgery, adni_design) and an arms selection.
#
# Usage (from a clean worktree at the run's pushed commit, submitted from the worktree's root; every path
# absolute; pre-create each profile's settings.json with the tier's arms so no two tasks race to write it, and
# mkdir logs before sbatch — SLURM opens --output before the script runs):
#   mkdir -p logs
#   # tier 1 (each tier on its own approval and cap)
#   ARMS=pls_cv,predictive_orthonormal,weighted
#   for P in committed conditions dispersion; do
#     PYTHONPATH=$PWD/src python scripts/candidate_simulation_check.py run --profile $P --arms $ARMS \
#       --out-dir /abs/run --shard 0/1 --limit 0
#   done
#   sbatch --array=0-399%64 --export=ALL,CSC_PROFILE=committed,CSC_ARMS=$ARMS,CSC_OUT=/abs/run,N_SHARDS=400 \
#     scripts/candidate_simulation_array.sbatch
#   # ... conditions and dispersion likewise; tier 2: surgery and adni_design; tier 3 appends sparse_block to the
#   # committed records with CSC_ARMS=sparse_block and the SAME N_SHARDS as tier 1 (a resume reads its own shard).
#   # after every task of a profile completes:
#   python scripts/candidate_simulation_check.py merge --profile committed --out-dir /abs/run
#   python scripts/candidate_simulation_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   CSC_PROFILE   committed | conditions | dispersion | surgery | adni_design.
#   CSC_ARMS      'all' or a comma list of pls_cv, predictive_orthonormal, weighted, sparse_block.
#   CSC_OUT       Absolute run root; shards land in $CSC_OUT/<profile>/records/.
#   N_SHARDS      Total shards = --array width.
# Optional:
#   CSC_ARGS      Extra run flags ("--limit 1" for timing). They must not change the settings between
#                 submissions: a resume into a directory with different settings is refused.
#   CSC_DRY_RUN   When set, print the run command instead of executing it.
#
# PYTHONPATH is set here to the submitting worktree's src/, so the run measures that revision's package.
# Each task runs its shard's units one at a time (one unit per CPU, BLAS pinned to one thread): a unit generates
# its dataset once, fits production once for pls_cv and predictive_orthonormal, and fits weighted and (on its
# three designs) sparse_block on their own, every arm at the unit's RRPP seed. Shards are round-robin over a
# design-major order and resumable by (unit, arm): resubmitting failed array ids skips complete keys. The `%64`
# throttle caps concurrent tasks (one node must not cold-start ~100 Python processes from /home1); change it
# mid-run with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${CSC_PROFILE:?CSC_PROFILE environment variable is required (committed, conditions, dispersion, surgery or adni_design).}"
: "${CSC_ARMS:?CSC_ARMS environment variable is required ('all' or a comma list of arms).}"
: "${CSC_OUT:?CSC_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

case "$CSC_PROFILE" in
    committed | conditions | dispersion | surgery | adni_design) ;;
    *)
        echo "CSC_PROFILE must be committed, conditions, dispersion, surgery or adni_design, got '$CSC_PROFILE'." >&2
        exit 2
        ;;
esac
case "$CSC_OUT" in
    /*) ;;
    *)
        echo "CSC_OUT must be an absolute path, got '$CSC_OUT'." >&2
        exit 2
        ;;
esac

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$REPO_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$REPO_DIR/src"

# shellcheck disable=SC2206
COMMAND=(python "$REPO_DIR/scripts/candidate_simulation_check.py" run --profile "$CSC_PROFILE" --arms "$CSC_ARMS"
         --out-dir "$CSC_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${CSC_ARGS:-})

if [ -n "${CSC_DRY_RUN:-}" ]; then
    printf 'PYTHONPATH=%q ' "$PYTHONPATH"
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

mkdir -p "$CSC_OUT" logs
if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
"${COMMAND[@]}"
