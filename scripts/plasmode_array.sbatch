#!/bin/bash
#SBATCH --job-name=motco-plasmode
#SBATCH --output=logs/motco-plasmode-%A_%a.out
#SBATCH --error=logs/motco-plasmode-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Plasmode check — paper-grade SLURM array template.
#
# Usage (from a worktree at the run's commit; every path absolute; the export copied and its sha256 verified
# against the settings; pre-create settings.json first so no two tasks race to write it):
#   mkdir -p logs
#   PYTHONPATH=$PWD/src python scripts/plasmode_check.py run --export-dir /abs/cohort_30 \
#     --fidelity /abs/calibration/fidelity.json --calibration /abs/calibration/calibration.json \
#     --out-dir /abs/run --index 0 --count 150 --limit 0
#   sbatch --array=0-149 --export=ALL,PYTHONPATH=$PWD/src,PM_EXPORT=/abs/cohort_30,\
#   PM_FIDELITY=/abs/calibration/fidelity.json,PM_CALIBRATION=/abs/calibration/calibration.json,\
#   PM_OUT=/abs/run,N_SHARDS=150 scripts/plasmode_array.sbatch
#   # after every task completes:
#   python scripts/plasmode_check.py merge --out-dir /abs/run
#   python scripts/plasmode_check.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   PM_EXPORT       The sead_sexdiff cohort_30 export directory.
#   PM_FIDELITY     fidelity.json of a passed gate on this export.
#   PM_CALIBRATION  calibration.json whose terminal iteration the settings adopt.
#   PM_OUT          Output directory; shards land in $PM_OUT/records/.
#   N_SHARDS        Total shards = --array width.
# Optional:
#   PM_ARGS         Extra run flags ("--limit 1" for timing). They must not change the settings between
#                   submissions: a resume into a directory with different settings is refused.
#
# Each unit runs one covariate-adjusted analysis with attribution (999 permutations, B = 200) and one unadjusted
# attribution on its frozen model. Each shard is resumable: resubmitting failed array ids skips recorded units.

set -euo pipefail

: "${PM_EXPORT:?PM_EXPORT environment variable is required.}"
: "${PM_FIDELITY:?PM_FIDELITY environment variable is required.}"
: "${PM_CALIBRATION:?PM_CALIBRATION environment variable is required.}"
: "${PM_OUT:?PM_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$PM_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/plasmode_check.py run --export-dir "$PM_EXPORT" --fidelity "$PM_FIDELITY" \
    --calibration "$PM_CALIBRATION" --out-dir "$PM_OUT" --index "${SLURM_ARRAY_TASK_ID}" --count "${N_SHARDS}" \
    ${PM_ARGS:-}
