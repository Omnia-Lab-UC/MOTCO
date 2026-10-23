#!/bin/bash
#SBATCH --job-name=motco-rank-floor
#SBATCH --output=logs/motco-rank-floor-%A_%a.out
#SBATCH --error=logs/motco-rank-floor-%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Rank-floor check — SLURM array template.
#
# Usage:
#   sbatch --array=0-49 \
#     --export=ALL,RANK_FLOOR_OUT=/abs/path/out,N_SHARDS=50,RANK_FLOOR_ARGS="--delta-methyl 0.5 --delta-expr 0.5" \
#     scripts/rank_floor_array.sbatch
#   # after every task completes:
#   python scripts/rank_floor_check.py merge --out-dir /abs/path/out
#
# Required environment variables (passed via --export):
#   RANK_FLOOR_OUT   Output directory; shards land in $RANK_FLOOR_OUT/records/.
#   N_SHARDS         Total shards = --array width.
# Optional:
#   RANK_FLOOR_ARGS  Extra run flags: the calibrated stage-signal strengths
#                    (e.g. "--delta-methyl 0.5 --delta-expr 0.5 --p-dmp 0.1"),
#                    or "--limit 1" for timing. They must not change the
#                    settings between submissions: a resume into a directory
#                    with different settings is refused.
#
# Each unit runs two analyses (floor_1, floor_2) at 999 permutations, about
# 2x a covariate Type I replicate. Each shard is resumable: resubmitting failed
# array ids skips recorded units.

set -euo pipefail

: "${RANK_FLOOR_OUT:?RANK_FLOOR_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$RANK_FLOOR_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/rank_floor_check.py run \
    --out-dir "$RANK_FLOOR_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${RANK_FLOOR_ARGS:-}
