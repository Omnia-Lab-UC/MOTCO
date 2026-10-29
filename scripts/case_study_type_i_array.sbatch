#!/bin/bash
#SBATCH --job-name=motco-case-study-type-i
#SBATCH --output=logs/motco-case-study-type-i-%A_%a.out
#SBATCH --error=logs/motco-case-study-type-i-%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Case-study configuration Type I check — SLURM array template.
#
# Usage:
#   sbatch --array=0-49 \
#     --export=ALL,CASE_STUDY_OUT=/abs/path/out,N_SHARDS=50,CASE_STUDY_ARGS="--delta 0.1 --s-z 3" \
#     scripts/case_study_type_i_array.sbatch
#   # after every task completes:
#   python scripts/case_study_type_i_check.py merge --out-dir /abs/path/out
#
# Required environment variables (passed via --export):
#   CASE_STUDY_OUT   Output directory; shards land in $CASE_STUDY_OUT/records/.
#   N_SHARDS         Total shards = --array width.
# Optional:
#   CASE_STUDY_ARGS  Extra run flags: the calibrated strengths (e.g.
#                    "--delta 0.1 --s-z 3"), or "--limit 1" for timing. They
#                    must not change the settings between submissions: a
#                    resume into a directory with different settings is refused.
#
# Each unit runs three analyses (adjusted, unadjusted, latent_stage) at 999
# permutations, about 3x a covariate Type I replicate. Each shard is
# resumable: resubmitting failed array ids skips recorded units, and a unit
# recorded under another shard count is skipped too.

set -euo pipefail

: "${CASE_STUDY_OUT:?CASE_STUDY_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$CASE_STUDY_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/case_study_type_i_check.py run \
    --out-dir "$CASE_STUDY_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${CASE_STUDY_ARGS:-}
