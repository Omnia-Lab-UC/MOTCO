#!/bin/bash
#SBATCH --job-name=motco-psc
#SBATCH --output=logs/motco-psc-%A_%a.out
#SBATCH --error=logs/motco-psc-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Predictive-subspace pilot (pls_cv, predictive, predictive_orthonormal, predictive_fitted, forced_k_minus_1,
# predictive_rank_plus_one) — SLURM array template for the three profiles.
#
# Usage (from a clean worktree at the run's pushed commit, which must carry MOTCO 0.24.0; every path absolute;
# pre-create each profile's settings.json first so no two tasks race to write it, and mkdir logs before
# sbatch — SLURM opens --output before the script runs):
#   mkdir -p logs
#   # simulation (7,800 units) and weak_size (1,000 units): the parents are read from results/
#   PYTHONPATH=$PWD/src python scripts/predictive_subspace_check.py run --profile simulation \
#     --out-dir /abs/run --shard 0/400 --limit 0
#   sbatch --array=0-399%64 --export=ALL,PSC_PROFILE=simulation,PSC_OUT=/abs/run,N_SHARDS=400 \
#     scripts/predictive_subspace_array.sbatch
#   PYTHONPATH=$PWD/src python scripts/predictive_subspace_check.py run --profile weak_size \
#     --out-dir /abs/run --shard 0/100 --limit 0
#   sbatch --array=0-99%64 --export=ALL,PSC_PROFILE=weak_size,PSC_OUT=/abs/run,N_SHARDS=100 \
#     scripts/predictive_subspace_array.sbatch
#   # plasmode (800 units): the export copied to the cluster; its file hashes must equal the plasmode's settings'
#   PYTHONPATH=$PWD/src python scripts/predictive_subspace_check.py run --profile plasmode \
#     --export-dir /abs/cohort_30 --out-dir /abs/run --shard 0/200 --limit 0
#   sbatch --array=0-199%64 --export=ALL,PSC_PROFILE=plasmode,PSC_OUT=/abs/run,PSC_EXPORT=/abs/cohort_30,\
#   N_SHARDS=200 scripts/predictive_subspace_array.sbatch
#   # after every task of the three arrays completes:
#   python scripts/predictive_subspace_check.py merge --profile simulation --out-dir /abs/run
#   python scripts/predictive_subspace_check.py merge --profile weak_size --out-dir /abs/run
#   python scripts/predictive_subspace_check.py merge --profile plasmode --out-dir /abs/run
#   python scripts/predictive_subspace_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   PSC_PROFILE   simulation | weak_size | plasmode.
#   PSC_OUT       Run root (absolute); shards land in $PSC_OUT/<profile>/records/.
#   N_SHARDS      Total shards = --array width.
# Plasmode only:
#   PSC_EXPORT    The sead_sexdiff cohort_30 export directory (absolute).
# Optional:
#   PSC_REPO      The worktree to run from (absolute; default the submission directory). Its src/ goes first on
#                 PYTHONPATH, so the run measures that worktree's package, not an installed one.
#   PSC_ARGS      Extra run flags ("--limit 1" for timing; "--reps 400" for the plasmode anchor as an approved
#                 lever). They must not change the settings between submissions: a resume into a directory with
#                 different settings is refused.
#   PSC_DRY_RUN   When set, print the run command (and the PYTHONPATH it would use) instead of executing it.
#
# Each task runs its shard's units one at a time (one CPU, BLAS pinned to one thread): a unit builds its dataset
# once, fits production once and evaluates every arm serially at the unit's RRPP seed. Shards are round-robin
# and resumable: resubmitting failed array ids skips complete units. The `%64` throttle caps concurrent tasks;
# raise it mid-run with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${PSC_PROFILE:?PSC_PROFILE environment variable is required (simulation, weak_size or plasmode).}"
: "${PSC_OUT:?PSC_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

REPO="${PSC_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"

for path_var in PSC_OUT REPO; do
    case "${!path_var}" in
        /*) ;;
        *)
            echo "$path_var must be an absolute path, got '${!path_var}'." >&2
            exit 2
            ;;
    esac
done

PROFILE_ARGS=()
case "$PSC_PROFILE" in
    simulation|weak_size) ;;
    plasmode)
        : "${PSC_EXPORT:?PSC_EXPORT environment variable is required for the plasmode profile.}"
        case "$PSC_EXPORT" in
            /*) ;;
            *)
                echo "PSC_EXPORT must be an absolute path, got '$PSC_EXPORT'." >&2
                exit 2
                ;;
        esac
        PROFILE_ARGS=(--export-dir "$PSC_EXPORT")
        ;;
    *)
        echo "PSC_PROFILE must be simulation, weak_size or plasmode, got '$PSC_PROFILE'." >&2
        exit 2
        ;;
esac

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

# shellcheck disable=SC2206
COMMAND=(python "$REPO/scripts/predictive_subspace_check.py" run --profile "$PSC_PROFILE" "${PROFILE_ARGS[@]}"
         --out-dir "$PSC_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${PSC_ARGS:-})

if [ -n "${PSC_DRY_RUN:-}" ]; then
    printf 'PYTHONPATH=%q ' "$PYTHONPATH"
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

cd "$REPO"
mkdir -p "$PSC_OUT" logs
if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
"${COMMAND[@]}"
