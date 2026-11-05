#!/bin/bash
#SBATCH --job-name=motco-lsc
#SBATCH --output=logs/motco-lsc-%A_%a.out
#SBATCH --error=logs/motco-lsc-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Latent-space comparison — SLURM array template.
#
# Usage (from a worktree at the run's commit; every path absolute; pre-create settings.json first so no two
# tasks race to write it, and mkdir logs before sbatch — SLURM opens --output before the script runs):
#   mkdir -p logs
#   PYTHONPATH=$PWD/src python scripts/latent_space_comparison.py run \
#     --calibration /abs/calibration/calibration.json --out-dir /abs/run --shard 0/200 --limit 0
#   sbatch --array=0-199%64 --export=ALL,PYTHONPATH=$PWD/src,LSC_CALIBRATION=/abs/calibration/calibration.json,\
#   LSC_OUT=/abs/run,N_SHARDS=200 scripts/latent_space_comparison_array.sbatch
#   # after every task completes:
#   python scripts/latent_space_comparison.py merge --out-dir /abs/run
#   python scripts/latent_space_comparison.py decide --out-dir /abs/run
#
# Required environment variables (passed via --export):
#   LSC_CALIBRATION  calibration.json whose adopted δ the weak design uses.
#   LSC_OUT          Output directory; shards land in $LSC_OUT/records/.
#   N_SHARDS         Total shards = --array width.
# Optional:
#   LSC_ARGS         Extra run flags ("--limit 1" for timing). They must not change the settings between
#                    submissions: a resume into a directory with different settings is refused.
#
# Each unit generates one dataset, fits stage-PLS once and evaluates five arms (pls_cv, pca_matched, pca_fixed,
# pls_residual_pca, concat) at one RRPP seed. Shards are round-robin over the 4,200 units and resumable:
# resubmitting failed array ids skips complete units. The `%64` throttle caps concurrent tasks; raise it mid-run
# with `scontrol update JobId=<array> ArrayTaskThrottle=<n>`.

set -euo pipefail

: "${LSC_CALIBRATION:?LSC_CALIBRATION environment variable is required.}"
: "${LSC_OUT:?LSC_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$LSC_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/latent_space_comparison.py run \
    --calibration "$LSC_CALIBRATION" \
    --out-dir "$LSC_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${LSC_ARGS:-}
