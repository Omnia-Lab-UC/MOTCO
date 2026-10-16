#!/bin/bash
#SBATCH --job-name=motco-cov
#SBATCH --output=logs/motco-cov-%A_%a.out
#SBATCH --error=logs/motco-cov-%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Covariate-adjustment Type I check — SLURM array template.
#
# Usage:
#   sbatch --array=0-49 \
#     --export=ALL,COV_OUT=/abs/path/out,N_SHARDS=50 \
#     scripts/covariate_type_i_array.sbatch
#   # after every task completes:
#   python scripts/covariate_type_i_check.py merge --out-dir /abs/path/out
#
# Required environment variables (passed via --export):
#   COV_OUT    Output directory; shards land in $COV_OUT/records/.
#   N_SHARDS   Total shards = --array width.
# Optional:
#   COV_ARGS   Extra run flags (e.g. "--limit 1" for timing, or calibrated
#              "--s-z 1.0"). They must not change the settings between
#              submissions: a resume into a directory with different settings
#              is refused.
#
# Each shard is resumable: resubmitting failed array ids skips recorded replicates.

set -euo pipefail

: "${COV_OUT:?COV_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$COV_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/covariate_type_i_check.py run \
    --out-dir "$COV_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${COV_ARGS:-}
