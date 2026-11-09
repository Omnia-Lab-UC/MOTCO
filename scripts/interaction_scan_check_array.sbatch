#!/bin/bash
#SBATCH --job-name=motco-isc
#SBATCH --output=logs/motco-isc-%A_%a.out
#SBATCH --error=logs/motco-isc-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Per-feature interaction scan validation — SLURM array template for both profiles.
#
# Usage (from a clean worktree at the run's commit; every path absolute; pre-create each profile's
# settings.json first so no two tasks race to write it, and mkdir logs before sbatch — SLURM opens --output
# before the script runs):
#   mkdir -p logs
#   # simulation: the split-test run's and the comparison's committed files are read from the worktree's results/
#   PYTHONPATH=$PWD/src python scripts/interaction_scan_check.py run --profile simulation \
#     --out-dir /abs/run --shard 0/300 --limit 0
#   sbatch --array=0-299%64 --export=ALL,PYTHONPATH=$PWD/src,ISC_PROFILE=simulation,ISC_OUT=/abs/run,\
#   N_SHARDS=300 scripts/interaction_scan_check_array.sbatch
#   # plasmode: the fidelity gate first, on the export copied to the cluster and verified by sha256
#   PYTHONPATH=$PWD/src python scripts/plasmode_check.py fidelity --export-dir /abs/cohort_30 \
#     --out-dir /abs/run/plasmode-fidelity
#   PYTHONPATH=$PWD/src python scripts/interaction_scan_check.py run --profile plasmode \
#     --export-dir /abs/cohort_30 --fidelity /abs/run/plasmode-fidelity/fidelity.json --out-dir /abs/run \
#     --shard 0/150 --limit 0
#   sbatch --array=0-149%64 --export=ALL,PYTHONPATH=$PWD/src,ISC_PROFILE=plasmode,ISC_OUT=/abs/run,\
#   ISC_EXPORT=/abs/cohort_30,ISC_FIDELITY=/abs/run/plasmode-fidelity/fidelity.json,N_SHARDS=150 \
#   scripts/interaction_scan_check_array.sbatch
#   # after every task of both arrays completes:
#   python scripts/interaction_scan_check.py merge --profile simulation --out-dir /abs/run
#   python scripts/interaction_scan_check.py merge --profile plasmode --out-dir /abs/run
#   python scripts/interaction_scan_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   ISC_PROFILE   simulation | plasmode.
#   ISC_OUT       Run root; shards land in $ISC_OUT/<profile>/records/.
#   N_SHARDS      Total shards = --array width.
# Plasmode only:
#   ISC_EXPORT    The sead_sexdiff cohort_30 export directory (file hashes must equal the plasmode's settings).
#   ISC_FIDELITY  fidelity.json of a passed gate on this export.
# Optional:
#   ISC_ARGS      Extra run flags ("--limit 1" for timing). They must not change the settings between
#                 submissions: a resume into a directory with different settings is refused.
#   ISC_DRY_RUN   When set, print the run command instead of executing it.
#
# Each task runs its shard's units one at a time (one unit's analysis per CPU, BLAS pinned to one thread):
# a simulation unit is the split-test validation's unit (the comparison's dataset, one pls_cv fit and
# evaluation, the split test) plus the scan on the split test's inputs; a plasmode unit is one
# analyze_trajectories with the split test and the scan. Shards are round-robin and resumable: resubmitting
# failed array ids skips complete units. The `%64` throttle caps concurrent tasks; raise it mid-run with
# `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${ISC_PROFILE:?ISC_PROFILE environment variable is required (simulation or plasmode).}"
: "${ISC_OUT:?ISC_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

mkdir -p "$ISC_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

PROFILE_ARGS=()
case "$ISC_PROFILE" in
    simulation) ;;
    plasmode)
        : "${ISC_EXPORT:?ISC_EXPORT environment variable is required for the plasmode profile.}"
        : "${ISC_FIDELITY:?ISC_FIDELITY environment variable is required for the plasmode profile.}"
        PROFILE_ARGS=(--export-dir "$ISC_EXPORT" --fidelity "$ISC_FIDELITY")
        ;;
    *)
        echo "ISC_PROFILE must be simulation or plasmode, got '$ISC_PROFILE'." >&2
        exit 2
        ;;
esac

# shellcheck disable=SC2206
COMMAND=(python scripts/interaction_scan_check.py run --profile "$ISC_PROFILE" "${PROFILE_ARGS[@]}"
         --out-dir "$ISC_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${ISC_ARGS:-})

if [ -n "${ISC_DRY_RUN:-}" ]; then
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi
"${COMMAND[@]}"
