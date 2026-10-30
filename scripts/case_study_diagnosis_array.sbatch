#!/bin/bash
#SBATCH --job-name=motco-case-study-diagnosis
#SBATCH --output=logs/motco-case-study-diagnosis-%A_%a.out
#SBATCH --error=logs/motco-case-study-diagnosis-%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Case-study configuration diagnosis — SLURM array template (one array per parent run).
#
# Usage:
#   sbatch --array=0-49 \
#     --export=ALL,DIAGNOSIS_OUT=/abs/path/out/delta0.1,DIAGNOSIS_PARENT=results/case-study-type-i-2026-10-29/delta0.1,N_SHARDS=50 \
#     scripts/case_study_diagnosis_array.sbatch
#   # after every task completes:
#   python scripts/case_study_diagnosis.py merge --out-dir /abs/path/out/delta0.1
#   # after both parents are merged:
#   python scripts/case_study_diagnosis.py decide --run /abs/path/out/delta0.1 --run /abs/path/out/delta0.05 \
#     --out-dir /abs/path/out
#
# Required environment variables (passed via --export):
#   DIAGNOSIS_OUT     Output directory for this parent; shards land in $DIAGNOSIS_OUT/records/.
#   DIAGNOSIS_PARENT  Committed run directory of the case-study Type I check (refused unless its settings
#                     pass their hash and the current code rebuilds them).
#   N_SHARDS          Total shards = --array width.
# Optional:
#   DIAGNOSIS_ARGS    Extra run flags: "--reps N" / "--perms B" (a prefix or a smoke run), or "--limit 1"
#                     for timing. They must not change the settings between submissions: a resume into a
#                     directory with different settings, or against another parent, is refused.
#
# Each unit runs four analyses (the covariate-blind basis with the six covariates, and the stage latent
# adjustment at CV floor 2, CV floor 1 and forced rank 2) plus 12 RRPP-only arms at the parent's 999
# permutations. Each shard is resumable: resubmitting failed array ids skips recorded units, and a unit
# recorded under another shard count is skipped too.

set -euo pipefail

: "${DIAGNOSIS_OUT:?DIAGNOSIS_OUT environment variable is required.}"
: "${DIAGNOSIS_PARENT:?DIAGNOSIS_PARENT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$DIAGNOSIS_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/case_study_diagnosis.py run \
    --parent "$DIAGNOSIS_PARENT" \
    --out-dir "$DIAGNOSIS_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${DIAGNOSIS_ARGS:-}
