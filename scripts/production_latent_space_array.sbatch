#!/bin/bash
#SBATCH --job-name=motco-plsc
#SBATCH --output=logs/motco-plsc-%A_%a.out
#SBATCH --error=logs/motco-plsc-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Production latent-space selection (pls_cv, weighted, sparse_block, weighted_sparse_block) — SLURM array
# template for both profiles.
#
# Usage (from a clean worktree at the run's commit, which must carry MOTCO 0.22.0; every path absolute;
# pre-create each profile's settings.json first so no two tasks race to write it, and mkdir logs before
# sbatch — SLURM opens --output before the script runs):
#   mkdir -p logs
#   # simulation: the comparison's committed settings, calibration and rate table are read from results/
#   PYTHONPATH=$PWD/src python scripts/production_latent_space_check.py run --profile simulation \
#     --out-dir /abs/run --shard 0/300 --limit 0
#   sbatch --array=0-299%64 --export=ALL,PYTHONPATH=$PWD/src,PLS_PROFILE=simulation,PLS_OUT=/abs/run,\
#   N_SHARDS=300 scripts/production_latent_space_array.sbatch
#   # plasmode: the export copied to the cluster; its file hashes must equal the plasmode's settings'
#   PYTHONPATH=$PWD/src python scripts/production_latent_space_check.py run --profile plasmode \
#     --export-dir /abs/cohort_30 --out-dir /abs/run --shard 0/100 --limit 0
#   sbatch --array=0-99%64 --export=ALL,PYTHONPATH=$PWD/src,PLS_PROFILE=plasmode,PLS_OUT=/abs/run,\
#   PLS_EXPORT=/abs/cohort_30,N_SHARDS=100 scripts/production_latent_space_array.sbatch
#   # after every task of both arrays completes:
#   python scripts/production_latent_space_check.py merge --profile simulation --out-dir /abs/run
#   python scripts/production_latent_space_check.py merge --profile plasmode --out-dir /abs/run
#   python scripts/production_latent_space_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   PLS_PROFILE   simulation | plasmode.
#   PLS_OUT       Run root (absolute); shards land in $PLS_OUT/<profile>/records/.
#   N_SHARDS      Total shards = --array width.
# Plasmode only:
#   PLS_EXPORT    The sead_sexdiff cohort_30 export directory (absolute).
# Optional:
#   PLS_ARGS      Extra run flags ("--limit 1" for timing; "--ladder 1,0.2,0.05" or "--plant-reps 50" as the
#                 approved levers). They must not change the settings between submissions: a resume into a
#                 directory with different settings is refused.
#   PLS_DRY_RUN   When set, print the run command instead of executing it.
#
# Each task runs its shard's units one at a time (one CPU, BLAS pinned to one thread): a unit builds its
# dataset once (the comparison's generator, or the plasmode's permuted labels and the span plant) and fits
# and evaluates the four arms serially at the unit's RRPP seed. Shards are round-robin and resumable:
# resubmitting failed array ids skips complete units. The `%64` throttle caps concurrent tasks; raise it
# mid-run with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${PLS_PROFILE:?PLS_PROFILE environment variable is required (simulation or plasmode).}"
: "${PLS_OUT:?PLS_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

case "$PLS_OUT" in
    /*) ;;
    *)
        echo "PLS_OUT must be an absolute path, got '$PLS_OUT'." >&2
        exit 2
        ;;
esac

PROFILE_ARGS=()
case "$PLS_PROFILE" in
    simulation) ;;
    plasmode)
        : "${PLS_EXPORT:?PLS_EXPORT environment variable is required for the plasmode profile.}"
        case "$PLS_EXPORT" in
            /*) ;;
            *)
                echo "PLS_EXPORT must be an absolute path, got '$PLS_EXPORT'." >&2
                exit 2
                ;;
        esac
        PROFILE_ARGS=(--export-dir "$PLS_EXPORT")
        ;;
    *)
        echo "PLS_PROFILE must be simulation or plasmode, got '$PLS_PROFILE'." >&2
        exit 2
        ;;
esac

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# shellcheck disable=SC2206
COMMAND=(python scripts/production_latent_space_check.py run --profile "$PLS_PROFILE" "${PROFILE_ARGS[@]}"
         --out-dir "$PLS_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${PLS_ARGS:-})

if [ -n "${PLS_DRY_RUN:-}" ]; then
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

mkdir -p "$PLS_OUT" logs
if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
"${COMMAND[@]}"
