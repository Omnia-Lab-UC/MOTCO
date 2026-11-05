#!/bin/bash
#SBATCH --job-name=motco-leak
#SBATCH --output=logs/motco-leak-%A_%a.out
#SBATCH --error=logs/motco-leak-%A_%a.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Stage-span leak decomposition — SLURM array template (scripts/stage_span_leak_check.py).
#
# Usage (from a worktree at the run's commit; every path absolute; mkdir logs before sbatch — SLURM opens
# --output before the script runs — and pass PYTHONPATH so the worktree's src is imported):
#   mkdir -p logs
#   sbatch --array=0-99%50 --export=ALL,PYTHONPATH=$PWD/src,LEAK_CONFIG=$PWD/examples/trajectory_power_study/\
#   stage_span_paper_grade_span_geometry.json,LEAK_OUT=/abs/leak/planted,N_SHARDS=100,\
#   LEAK_ARGS="--span-scale planted --modes orientation shape --effects 0.1 0.25 0.5 --replicates 100 --permutations 199" \
#   scripts/stage_span_leak_check_array.sbatch
#   # after every task completes, per arm:
#   python scripts/stage_span_leak_check.py merge --out-dir /abs/leak/planted
#   python scripts/stage_span_leak_check.py gate --planted /abs/leak/planted --null /abs/leak/null \
#     --out /abs/leak/preflight_gate.json
#
# Required environment variables (passed via --export):
#   LEAK_CONFIG  The profile whose stage-span cells are decomposed.
#   LEAK_OUT     Output directory for shard_<i>.csv and settings.json.
#   N_SHARDS     Total shards = --array width.
#   LEAK_ARGS    The run flags (--span-scale, --columns, --modes, --effects, --replicates, --permutations). They
#                must not change between submissions: a resume into a directory with other settings is refused.
#
# Shards are round-robin over the units and resumable: resubmitting failed array ids skips complete units.

set -euo pipefail

: "${LEAK_CONFIG:?LEAK_CONFIG environment variable is required.}"
: "${LEAK_OUT:?LEAK_OUT environment variable is required.}"
: "${N_SHARDS:?N_SHARDS environment variable is required.}"
: "${LEAK_ARGS:?LEAK_ARGS environment variable is required.}"

mkdir -p "$LEAK_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi

# shellcheck disable=SC2086
python scripts/stage_span_leak_check.py run \
    --config "$LEAK_CONFIG" \
    --out-dir "$LEAK_OUT" \
    --shard "$SLURM_ARRAY_TASK_ID" \
    --n-shards "$N_SHARDS" \
    $LEAK_ARGS
