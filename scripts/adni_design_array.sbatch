#!/bin/bash
#SBATCH --job-name=motco-adni-design
#SBATCH --output=logs/motco-adni-design-%A_%a.out
#SBATCH --error=logs/motco-adni-design-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# ADNI-design operating characteristics — SLURM array template (one array for the paper-grade run).
#
# Usage:
#   # pre-create settings.json so concurrent shards never race on it
#   python scripts/adni_design_check.py run --shard 0/200 --limit 0 --out-dir /abs/path/out/run
#   sbatch --array=0-199 --export=ALL,ADNI_OUT=/abs/path/out/run,N_SHARDS=200 scripts/adni_design_array.sbatch
#   # after every task completes:
#   python scripts/adni_design_check.py merge --out-dir /abs/path/out/run
#   python scripts/adni_design_check.py decide --out-dir /abs/path/out/run
#
# Required environment variables (passed via --export):
#   ADNI_OUT      Output directory of the run; shards land in $ADNI_OUT/records/.
#   N_SHARDS      Total shards = --array width.
# Optional:
#   ADNI_ARGS     Extra run flags: the calibrated strengths ("--delta D --s-z Z --s-p P") when they differ from the
#                 script's defaults, "--reps N" (the approved replicate count), "--perms B" (smoke), or "--limit 1"
#                 for timing. They must not change the settings between submissions: a resume into a directory with
#                 different settings is refused.
#
# Each unit is one adjusted analysis at n = 1,501 and 999 permutations. Each shard is resumable: resubmitting
# failed or stalled array ids skips recorded units, and a unit recorded under another shard count is skipped too.

set -euo pipefail

: "${ADNI_OUT:?ADNI_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$ADNI_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/adni_design_check.py run --out-dir "$ADNI_OUT" --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" ${ADNI_ARGS:-}
