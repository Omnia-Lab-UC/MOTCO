#!/bin/bash
#SBATCH --job-name=motco-gain
#SBATCH --output=logs/motco-gain-%A_%a.out
#SBATCH --error=logs/motco-gain-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# ADNI gain and specificity — SLURM array template for the simulation profile only (stage-span plants at the
# ADNI-design simulation; four arms per unit on one fit). The plants profile reads ADNI's export, which is never
# copied between machines, so it has no array: it runs on the workstation where the export is.
#
# Usage (from a clean worktree at the run's pushed commit; every path absolute; mkdir logs before sbatch — SLURM
# opens --output before the script runs). Run the simulation pre-flight once first, so no two tasks race to write
# the settings and every task refuses to start without a passed pre-flight:
#   mkdir -p logs
#   PYTHONPATH=$PWD/src python scripts/adni_gain_specificity_check.py preflight --profile simulation --out-dir /abs/run
#   sbatch --array=0-99%50 --export=ALL,GAIN_PROFILE=simulation,GAIN_OUT=/abs/run,N_SHARDS=100 \
#     scripts/adni_gain_specificity_array.sbatch
#   # after every task completes:
#   python scripts/adni_gain_specificity_check.py merge --profile simulation --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   GAIN_PROFILE  simulation (plants is refused: the export is not copied between machines).
#   GAIN_OUT      Run root (absolute); shards land in $GAIN_OUT/simulation/records/.
#   N_SHARDS      Total shards = --array width.
# Optional:
#   GAIN_REPO     The worktree to run from (absolute; default the submission directory). Its src/ goes first on
#                 PYTHONPATH, so the run measures that worktree's package, not an installed one.
#   GAIN_ARGS     Extra run flags ("--limit 1" for timing; "--control-reps 2 --plant-reps 2" for the timing shard).
#                 They must not change the settings between submissions: a resume into a directory with different
#                 settings is refused.
#   GAIN_DRY_RUN  When set, print the run command (and the PYTHONPATH it would use) instead of executing it.
#
# Each task runs its shard's units one at a time (one CPU, BLAS pinned to one thread). Shards are round-robin and
# resumable: resubmitting failed array ids skips complete units. The `%N` throttle caps concurrent tasks; raise it
# mid-run with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${GAIN_PROFILE:?GAIN_PROFILE environment variable is required.}"
: "${GAIN_OUT:?GAIN_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array).}"

case "$GAIN_PROFILE" in
    simulation) ;;
    plants)
        echo "GAIN_PROFILE=plants is refused: the plants profile reads ADNI's export, which is never copied between" \
             "machines; run it on the workstation where the export is." >&2
        exit 2
        ;;
    *)
        echo "GAIN_PROFILE must be simulation, got '$GAIN_PROFILE'." >&2
        exit 2
        ;;
esac

REPO="${GAIN_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"

for path_var in GAIN_OUT REPO; do
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
COMMAND=(python "$REPO/scripts/adni_gain_specificity_check.py" run --profile "$GAIN_PROFILE" --out-dir "$GAIN_OUT"
         --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${GAIN_ARGS:-})

if [ -n "${GAIN_DRY_RUN:-}" ]; then
    printf 'PYTHONPATH=%q ' "$PYTHONPATH"
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
    exit 0
fi

cd "$REPO"
mkdir -p "$GAIN_OUT" logs
if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
"${COMMAND[@]}"
