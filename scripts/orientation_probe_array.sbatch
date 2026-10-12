#!/bin/bash
#SBATCH --job-name=motco-omp
#SBATCH --output=logs/motco-omp-%A_%a.out
#SBATCH --error=logs/motco-omp-%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Orientation measurement probe — SLURM array template.
#
# Usage:
#   sbatch --array=0-49 \
#     --export=ALL,PROBE_OUT=/abs/path/out,N_SHARDS=50 \
#     scripts/orientation_probe_array.sbatch
#   # after every task completes:
#   python scripts/orientation_measurement_probe.py --merge --out-dir /abs/path/out
#
# Required environment variables (passed via --export):
#   PROBE_OUT  Output directory; shards land in $PROBE_OUT/records/.
#   N_SHARDS   Total shards = --array width.
# Optional:
#   PROBE_ARGS Extra probe flags (e.g. "--limit 1" for timing). They must
#              not change the probe settings between submissions: a resume
#              into a directory with different settings is refused.
#
# Each shard is resumable: resubmitting failed array ids skips completed units.

set -euo pipefail

: "${PROBE_OUT:?PROBE_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"

mkdir -p "$PROBE_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/orientation_measurement_probe.py \
    --out-dir "$PROBE_OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}/${N_SHARDS}" \
    ${PROBE_ARGS:-}
