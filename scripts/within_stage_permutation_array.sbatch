#!/bin/bash
#SBATCH --job-name=motco-wsp
#SBATCH --output=logs/motco-wsp-%A_%a.out
#SBATCH --error=logs/motco-wsp-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Within-stage permutation pilot (eight arms per fitted unit: the scores and the three predictive coordinate
# systems, each under the unrestricted and the within-stage draws; two arms per dispersion unit) — SLURM array
# template, one profile per array.
#
# Usage (from a clean worktree at the run's pushed commit; every path absolute; mkdir logs before sbatch — SLURM
# opens --output before the script runs). For the plasmode profile, copy the export to the cluster (its file
# hashes equal to the plasmode's settings') and run the pre-flight once first, so no two tasks race to write the
# settings and every task refuses to start without a passed pre-flight:
#   mkdir -p logs
#   sbatch --array=0-99%50 --export=ALL,WSP_PROFILE=simulation,WSP_OUT=/abs/run,N_SHARDS=100 \
#     scripts/within_stage_permutation_array.sbatch
#   sbatch --array=0-19%20 --export=ALL,WSP_PROFILE=weak_size,WSP_OUT=/abs/run,N_SHARDS=20 \
#     scripts/within_stage_permutation_array.sbatch
#   PYTHONPATH=$PWD/src python scripts/within_stage_permutation_check.py preflight \
#     --export-dir /abs/cohort_30 --out-dir /abs/run
#   sbatch --array=0-9%10 --export=ALL,WSP_PROFILE=plasmode,WSP_OUT=/abs/run,WSP_EXPORT=/abs/cohort_30,N_SHARDS=10 \
#     scripts/within_stage_permutation_array.sbatch
#   # after every task of a profile completes:
#   python scripts/within_stage_permutation_check.py merge --profile <profile> --out-dir /abs/run
#   # after every profile is merged:
#   python scripts/within_stage_permutation_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   WSP_PROFILE   simulation | weak_size | plasmode | dispersion.
#   WSP_OUT       Run root (absolute); shards land in $WSP_OUT/<profile>/records/.
#   N_SHARDS      Total shards = --array width.
# Required for the plasmode profile:
#   WSP_EXPORT    The sead_sexdiff cohort_30 export directory (absolute).
# Optional:
#   WSP_REPO      The worktree to run from (absolute; default the submission directory). Its src/ goes first on
#                 PYTHONPATH, so the run measures that worktree's package, not an installed one.
#   WSP_ARGS      Extra run flags ("--limit 1" for timing). They must not change the settings between submissions:
#                 a resume into a directory with different settings is refused.
#   WSP_DRY_RUN   When set, print the run command (and the PYTHONPATH it would use) instead of executing it.
#
# Each task runs its shard's units one at a time (one CPU, BLAS pinned to one thread). Shards are round-robin and
# resumable: resubmitting failed array ids skips complete units. The `%N` throttle caps concurrent tasks; raise it
# mid-run with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${WSP_PROFILE:?WSP_PROFILE environment variable is required.}"
: "${WSP_OUT:?WSP_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

case "$WSP_PROFILE" in
    simulation|weak_size|dispersion) EXPORT_ARGS=() ;;
    plasmode)
        : "${WSP_EXPORT:?WSP_EXPORT environment variable is required for the plasmode profile.}"
        EXPORT_ARGS=(--export-dir "$WSP_EXPORT")
        ;;
    *)
        echo "WSP_PROFILE must be simulation, weak_size, plasmode or dispersion, got '$WSP_PROFILE'." >&2
        exit 2
        ;;
esac

REPO="${WSP_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"

for path_var in WSP_OUT REPO; do
    case "${!path_var}" in
        /*) ;;
        *)
            echo "$path_var must be an absolute path, got '${!path_var}'." >&2
            exit 2
            ;;
    esac
done
if [ "$WSP_PROFILE" = plasmode ]; then
    case "$WSP_EXPORT" in
        /*) ;;
        *)
            echo "WSP_EXPORT must be an absolute path, got '$WSP_EXPORT'." >&2
            exit 2
            ;;
    esac
fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

# shellcheck disable=SC2206
COMMAND=(python "$REPO/scripts/within_stage_permutation_check.py" run --profile "$WSP_PROFILE"
         ${EXPORT_ARGS[@]+"${EXPORT_ARGS[@]}"} --out-dir "$WSP_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}"
         ${WSP_ARGS:-})

if [ -n "${WSP_DRY_RUN:-}" ]; then
    printf 'PYTHONPATH=%q ' "$PYTHONPATH"
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

cd "$REPO"
mkdir -p "$WSP_OUT" logs
if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
"${COMMAND[@]}"
