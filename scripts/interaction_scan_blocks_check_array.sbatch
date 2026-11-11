#!/bin/bash
#SBATCH --job-name=motco-isb
#SBATCH --output=logs/motco-isb-%A_%a.out
#SBATCH --error=logs/motco-isb-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Per-feature interaction scan by block — SLURM array template for the three profiles.
#
# Usage (from a clean worktree at the run's commit; every path absolute; pre-create each profile's
# settings.json first so no two tasks race to write it, and mkdir logs before sbatch — SLURM opens --output
# before the script runs):
#   mkdir -p logs
#   # plasmode profiles: the export copied to the cluster and verified by sha256 (no fidelity gate: the scan
#   # never reads the latent space, and the export is bound by its hashes)
#   PYTHONPATH=$PWD/src python scripts/interaction_scan_blocks_check.py run --profile plasmode \
#     --export-dir /abs/cohort_30 --out-dir /abs/run --shard 0/50 --limit 0
#   sbatch --array=0-49%50 --export=ALL,PYTHONPATH=$PWD/src,ISB_PROFILE=plasmode,ISB_OUT=/abs/run,\
#   ISB_EXPORT=/abs/cohort_30,N_SHARDS=50 scripts/interaction_scan_blocks_check_array.sbatch
#   # the same for --profile plasmode_unequal (ISB_PROFILE=plasmode_unequal)
#   # simulation (normally on the workstation, but the template takes it too): no export
#   sbatch --array=0-99%64 --export=ALL,PYTHONPATH=$PWD/src,ISB_PROFILE=simulation,ISB_OUT=/abs/run,\
#   N_SHARDS=100 scripts/interaction_scan_blocks_check_array.sbatch
#   # after every task completes:
#   python scripts/interaction_scan_blocks_check.py merge --profile <profile> --out-dir /abs/run   # each profile
#   python scripts/interaction_scan_blocks_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   ISB_PROFILE   simulation | plasmode | plasmode_unequal.
#   ISB_OUT       Run root; shards land in $ISB_OUT/<profile>/records/.
#   N_SHARDS      Total shards = --array width.
# Plasmode profiles only:
#   ISB_EXPORT    The sead_sexdiff cohort_30 export directory (file hashes must equal the plasmode's settings).
# Optional:
#   ISB_ARGS      Extra run flags ("--limit 1" for timing). They must not change the settings between
#                 submissions: a resume into a directory with different settings is refused.
#   ISB_DRY_RUN   When set, print the run command instead of executing it.
#
# Each task runs its shard's units one at a time (BLAS pinned to one thread): a unit builds the scan's inputs
# without fitting a latent space and runs the joint and by-block scans on them. Shards are round-robin and
# resumable: resubmitting failed array ids skips complete units. The `%N` throttle caps concurrent tasks; raise
# it mid-run with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${ISB_PROFILE:?ISB_PROFILE environment variable is required (simulation, plasmode or plasmode_unequal).}"
: "${ISB_OUT:?ISB_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

mkdir -p "$ISB_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

PROFILE_ARGS=()
case "$ISB_PROFILE" in
    simulation) ;;
    plasmode|plasmode_unequal)
        : "${ISB_EXPORT:?ISB_EXPORT environment variable is required for the plasmode profiles.}"
        PROFILE_ARGS=(--export-dir "$ISB_EXPORT")
        ;;
    *)
        echo "ISB_PROFILE must be simulation, plasmode or plasmode_unequal, got '$ISB_PROFILE'." >&2
        exit 2
        ;;
esac

# shellcheck disable=SC2206
COMMAND=(python scripts/interaction_scan_blocks_check.py run --profile "$ISB_PROFILE" "${PROFILE_ARGS[@]}"
         --out-dir "$ISB_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${ISB_ARGS:-})

if [ -n "${ISB_DRY_RUN:-}" ]; then
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi
"${COMMAND[@]}"
