#!/bin/bash
#SBATCH --job-name=motco-corrected-permutation
#SBATCH --output=logs/motco-corrected-permutation-%A_%a.out
#SBATCH --error=logs/motco-corrected-permutation-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Corrected-permutation check — SLURM array template (one array per run).
#
# Usage (the four Type I runs and the power run):
#   for table in astrocyte microglia; do for delta in 0.1 0.05; do
#     sbatch --array=0-49 \
#       --export=ALL,CP_OUT=/abs/path/out/type-i-${table}-delta${delta},CP_PROFILE=type_i,CP_TABLE=${table},CP_DELTA=${delta},N_SHARDS=50 \
#       scripts/corrected_permutation_array.sbatch
#   done; done
#   sbatch --array=0-139 \
#     --export=ALL,CP_OUT=/abs/path/out/power,CP_PROFILE=power,N_SHARDS=140 \
#     scripts/corrected_permutation_array.sbatch
#   # after every task of a run completes:
#   python scripts/corrected_permutation_check.py merge --out-dir /abs/path/out/<run>
#   # after the five runs are merged:
#   python scripts/corrected_permutation_check.py decide --runs /abs/path/out/type-i-* /abs/path/out/power \
#     --out-dir /abs/path/out
#
# Required environment variables (passed via --export):
#   CP_OUT        Output directory of this run; shards land in $CP_OUT/records/.
#   CP_PROFILE    type_i or power.
#   N_SHARDS      Total shards = --array width.
# Optional:
#   CP_TABLE      astrocyte (default) or microglia (type_i only; power runs the astrocyte table).
#   CP_DELTA      Per-omic stage effect (type_i: 0.1 or 0.05; power: 2.0, the default).
#   CP_ARGS       Extra run flags: "--reps N" / "--perms B" (smoke), or "--limit 1" for timing. They must not
#                 change the settings between submissions: a resume into a directory with different settings
#                 is refused.
#
# Each unit runs three analyses (no correction, df, leverage) at 999 permutations on the same data, CV
# settings and RRPP seed. Each shard is resumable: resubmitting failed array ids skips recorded units, and a
# unit recorded under another shard count is skipped too.

set -euo pipefail

: "${CP_OUT:?CP_OUT environment variable is required.}"
: "${CP_PROFILE:?CP_PROFILE environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$CP_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

run_args=(--profile "$CP_PROFILE" --out-dir "$CP_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}")
if [ -n "${CP_TABLE:-}" ]; then
    run_args+=(--table "$CP_TABLE")
fi
if [ -n "${CP_DELTA:-}" ]; then
    run_args+=(--delta "$CP_DELTA")
fi

# shellcheck disable=SC2086
python scripts/corrected_permutation_check.py run "${run_args[@]}" ${CP_ARGS:-}
