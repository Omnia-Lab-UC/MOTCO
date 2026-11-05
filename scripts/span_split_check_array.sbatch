#!/bin/bash
#SBATCH --job-name=motco-ssc
#SBATCH --output=logs/motco-ssc-%A_%a.out
#SBATCH --error=logs/motco-ssc-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Split interaction test validation — SLURM array template for both profiles.
#
# Usage (from a clean worktree at the run's commit; every path absolute; pre-create each profile's
# settings.json first so no two tasks race to write it, and mkdir logs before sbatch — SLURM opens --output
# before the script runs):
#   mkdir -p logs
#   # simulation: the comparison's calibration and rate table are read from the worktree's results/
#   PYTHONPATH=$PWD/src python scripts/span_split_check.py run --profile simulation \
#     --out-dir /abs/run --shard 0/300 --limit 0
#   sbatch --array=0-299%64 --export=ALL,PYTHONPATH=$PWD/src,SSC_PROFILE=simulation,SSC_OUT=/abs/run,\
#   N_SHARDS=300 scripts/span_split_check_array.sbatch
#   # plasmode: the fidelity gate first, on the export copied to the cluster and verified by sha256
#   PYTHONPATH=$PWD/src python scripts/plasmode_check.py fidelity --export-dir /abs/cohort_30 \
#     --out-dir /abs/run/plasmode-fidelity
#   PYTHONPATH=$PWD/src python scripts/span_split_check.py run --profile plasmode --export-dir /abs/cohort_30 \
#     --fidelity /abs/run/plasmode-fidelity/fidelity.json --out-dir /abs/run --shard 0/150 --limit 0
#   sbatch --array=0-149%64 --export=ALL,PYTHONPATH=$PWD/src,SSC_PROFILE=plasmode,SSC_OUT=/abs/run,\
#   SSC_EXPORT=/abs/cohort_30,SSC_FIDELITY=/abs/run/plasmode-fidelity/fidelity.json,N_SHARDS=150 \
#   scripts/span_split_check_array.sbatch
#   # after every task of both arrays completes:
#   python scripts/span_split_check.py merge --profile simulation --out-dir /abs/run
#   python scripts/span_split_check.py merge --profile plasmode --out-dir /abs/run
#   python scripts/span_split_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   SSC_PROFILE   simulation | plasmode.
#   SSC_OUT       Run root; shards land in $SSC_OUT/<profile>/records/.
#   N_SHARDS      Total shards = --array width.
# Plasmode only:
#   SSC_EXPORT    The sead_sexdiff cohort_30 export directory (file hashes must equal the plasmode's settings).
#   SSC_FIDELITY  fidelity.json of a passed gate on this export.
# Optional:
#   SSC_ARGS      Extra run flags ("--limit 1" for timing). They must not change the settings between
#                 submissions: a resume into a directory with different settings is refused.
#
# Each task runs its shard's units one at a time (one unit's analysis per CPU, BLAS pinned to one thread):
# a simulation unit generates the comparison's dataset, fits pls_cv once, evaluates it and runs the split test;
# a plasmode unit is one analyze_trajectories with the split test. Shards are round-robin and resumable:
# resubmitting failed array ids skips complete units. The `%64` throttle caps concurrent tasks; raise it mid-run
# with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${SSC_PROFILE:?SSC_PROFILE environment variable is required (simulation or plasmode).}"
: "${SSC_OUT:?SSC_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$SSC_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

PROFILE_ARGS=()
case "$SSC_PROFILE" in
    simulation) ;;
    plasmode)
        : "${SSC_EXPORT:?SSC_EXPORT environment variable is required for the plasmode profile.}"
        : "${SSC_FIDELITY:?SSC_FIDELITY environment variable is required for the plasmode profile.}"
        PROFILE_ARGS=(--export-dir "$SSC_EXPORT" --fidelity "$SSC_FIDELITY")
        ;;
    *)
        echo "SSC_PROFILE must be simulation or plasmode, got '$SSC_PROFILE'." >&2
        exit 2
        ;;
esac

# shellcheck disable=SC2086
python scripts/span_split_check.py run --profile "$SSC_PROFILE" "${PROFILE_ARGS[@]}" \
    --out-dir "$SSC_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${SSC_ARGS:-}
