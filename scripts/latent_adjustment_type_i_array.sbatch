#!/bin/bash
#SBATCH --job-name=motco-latent
#SBATCH --output=logs/motco-latent-%A_%a.out
#SBATCH --error=logs/motco-latent-%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Latent-adjustment Type I check — SLURM array template.
#
# Usage:
#   sbatch --array=0-49 \
#     --export=ALL,LATENT_OUT=/abs/path/out,N_SHARDS=50 \
#     scripts/latent_adjustment_type_i_array.sbatch
#   # after every task completes:
#   python scripts/latent_adjustment_type_i_check.py merge --out-dir /abs/path/out
#
# Required environment variables (passed via --export):
#   LATENT_OUT   Output directory; shards land in $LATENT_OUT/records/.
#   N_SHARDS     Total shards = --array width.
# Optional:
#   LATENT_ARGS  Extra run flags (e.g. "--limit 1" for timing, or the calibrated
#                "--s-z 20"). They must not change the settings between
#                submissions: a resume into a directory with different settings
#                is refused.
#
# Each replicate runs three analyses (baseline, stage, group_stage), about 1.5x
# a covariate Type I replicate. Each shard is resumable: resubmitting failed
# array ids skips recorded replicates.

set -euo pipefail

: "${LATENT_OUT:?LATENT_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$LATENT_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/latent_adjustment_type_i_check.py run \
    --out-dir "$LATENT_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${LATENT_ARGS:-}
