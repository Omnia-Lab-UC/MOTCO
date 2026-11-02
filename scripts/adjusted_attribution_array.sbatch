#!/bin/bash
#SBATCH --job-name=motco-adjusted-attribution
#SBATCH --output=logs/motco-adjusted-attribution-%A_%a.out
#SBATCH --error=logs/motco-adjusted-attribution-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Adjusted-attribution check — SLURM array template.
#
# Usage (the paper-grade run; pre-create settings.json first so no two tasks race to write it):
#   python scripts/adjusted_attribution_check.py run --s-z 10 \
#     --calibration results/adjusted-attribution-2026-11-02/calibration/calibration.json \
#     --out-dir /abs/path/out/run --shard 0/50 --limit 0
#   sbatch --array=0-49 \
#     --export=ALL,AA_OUT=/abs/path/out/run,N_SHARDS=50,AA_CALIBRATION=results/adjusted-attribution-2026-11-02/calibration/calibration.json \
#     scripts/adjusted_attribution_array.sbatch
#   # after every task completes:
#   python scripts/adjusted_attribution_check.py merge --out-dir /abs/path/out/run
#   python scripts/adjusted_attribution_check.py decide --out-dir /abs/path/out/run
#
# Required environment variables (passed via --export):
#   AA_OUT           Output directory; shards land in $AA_OUT/records/.
#   N_SHARDS         Total shards = --array width.
#   AA_CALIBRATION   The calibration.json whose adopted iteration the settings name.
# Optional:
#   AA_S_Z           Injection strength (default 10, the adopted calibration iteration's).
#   AA_ARGS          Extra run flags: "--reps N" / "--perms P" / "--bootstrap B" (smoke), or "--limit 1" for timing.
#                    They must not change the settings between submissions: a resume into a directory with
#                    different settings is refused.
#
# Each unit runs one covariate-adjusted analysis with attribution (999 permutations, B = 200) and two more
# attributions on its frozen model (oracle and unadjusted). Each shard is resumable: resubmitting failed array
# ids skips recorded units, and a unit recorded under another shard count is skipped too.

set -euo pipefail

: "${AA_OUT:?AA_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${AA_CALIBRATION:?AA_CALIBRATION environment variable is required.}"

mkdir -p "$AA_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/adjusted_attribution_check.py run --s-z "${AA_S_Z:-10}" --calibration "$AA_CALIBRATION" \
    --out-dir "$AA_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${AA_ARGS:-}
