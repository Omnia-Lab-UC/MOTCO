#!/bin/bash
#SBATCH --job-name=motco-plasmode-calibrate
#SBATCH --output=logs/motco-plasmode-calibrate-%j.out
#SBATCH --error=logs/motco-plasmode-calibrate-%j.err
#SBATCH --partition=512x1024
#SBATCH --exclude=n11
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G

# Plasmode check — fidelity gate and strength calibration as one multi-CPU job.
#
# Usage (from a worktree at the run's commit, the export copied and its sha256 verified):
#   mkdir -p logs
#   sbatch --export=ALL,PYTHONPATH=$PWD/src,PM_EXPORT=/abs/path/cohort_30,PM_OUT=/abs/path/out/calibration \
#     scripts/plasmode_calibrate.sbatch
#
# Required environment variables (passed via --export):
#   PM_EXPORT   The sead_sexdiff cohort_30 export directory (read in place).
#   PM_OUT      Output directory: fidelity.json, then calibration.json (one entry per ladder step).
# Optional:
#   PM_WORKERS  Units evaluated in parallel (default: the allocated CPUs; 40 = one round per ladder step).
#
# The fidelity gate runs first (skipped when PM_OUT already holds a passed fidelity.json) and stops the job on a
# mismatch. Calibration then climbs e in {0.5, 1.0, 1.5, 2.0, 3.0}, 40 units per step at 99 permutations with
# attribution off, and resumes after the last recorded step when resubmitted.

set -euo pipefail

: "${PM_EXPORT:?PM_EXPORT environment variable is required.}"
: "${PM_OUT:?PM_OUT environment variable is required.}"

mkdir -p "$PM_OUT" logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -f .venv/bin/activate ]; then
    # shellcheck source=/dev/null
    source .venv/bin/activate
fi
python -c "import motco; print('motco', motco.__version__, motco.__file__)"

if [ ! -f "$PM_OUT/fidelity.json" ]; then
    python scripts/plasmode_check.py fidelity --export-dir "$PM_EXPORT" --out-dir "$PM_OUT" --n-jobs 5
fi
python scripts/plasmode_check.py calibrate --export-dir "$PM_EXPORT" --fidelity "$PM_OUT/fidelity.json" \
    --out-dir "$PM_OUT" --workers "${PM_WORKERS:-${SLURM_CPUS_PER_TASK:-1}}"
