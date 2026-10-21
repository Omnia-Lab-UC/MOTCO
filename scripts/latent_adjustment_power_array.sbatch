#!/bin/bash
#SBATCH --job-name=motco-latent-power
#SBATCH --output=logs/motco-latent-power-%A_%a.out
#SBATCH --error=logs/motco-latent-power-%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Latent-adjustment power comparison — SLURM array template.
#
# Usage:
#   sbatch --array=0-199 \
#     --export=ALL,POWER_OUT=/abs/path/out,N_SHARDS=200 \
#     scripts/latent_adjustment_power_array.sbatch
#   # after every task completes:
#   python scripts/latent_adjustment_power_comparison.py merge --out-dir /abs/path/out
#
# Required environment variables (passed via --export):
#   POWER_OUT   Output directory; shards land in $POWER_OUT/records/.
#   N_SHARDS    Total shards = --array width.
# Optional:
#   POWER_ARGS  Extra run flags (e.g. "--limit 1" for timing, or the settled
#               "--magnitude-effects 0.02,0.05,0.1,0.25"). They must not change
#               the settings between submissions: a resume into a directory with
#               different settings is refused.
#
# Units are (column, cell, replicate); shard i takes every N-th unit, so shards
# balance across cells. Each unit runs three analyses (baseline, stage,
# group_stage). Each shard is resumable: resubmitting failed array ids skips
# recorded units.

set -euo pipefail

: "${POWER_OUT:?POWER_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$POWER_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/latent_adjustment_power_comparison.py run \
    --out-dir "$POWER_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${POWER_ARGS:-}
