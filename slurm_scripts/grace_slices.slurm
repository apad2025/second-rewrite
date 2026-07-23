#!/bin/bash
#############################################################################
# grace_slices.slurm         Stage 2/3 of the slice-wise separation pipeline
#
# Separates one slice per array task. Each task is a plain single-process
# MATLAB, so no parallel pool and no MDCS licenses are involved - which is
# what the parfor version could never get past on Grace.
#
# Submit using the command:   sbatch grace_slices.slurm
# Check progress:             squeue -u $USER   /   ls FWbigc_zip_slices/
# Rerun a failed slice:       sbatch --array=37 grace_slices.slurm
#
# Timing note: slice 12 took ~712 s when run on its own, so 2 h leaves plenty
# of margin. Adjust --array to the slice count (50 for the SUSHI dates).
# The %10 throttle keeps concurrent MATLAB license checkouts reasonable;
# raise it if `lmstat -f MATLAB` shows spare seats.
#############################################################################
#SBATCH --job-name=fwsep_slice
#SBATCH --array=1-50%10
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=fwsep_slice_%A_%a.out
#SBATCH --error=fwsep_slice_%A_%a.err
### ---- EDIT THESE ---------------------------------------------------------
#SBATCH --account=ACCOUNT_XXXXXX
#SBATCH --partition=medium
#############################################################################

set -euo pipefail

MATLAB_MODULE="Matlab/R2024a"
CODE_DIR="/scratch/user/apad/GRMD/Images"

module purge
module load "$MATLAB_MODULE"

# One preferences dir per array task: a shared dir races once tasks overlap
export MATLAB_PREFDIR="$SCRATCH/mathworks_cache/prefdir_$SLURM_ARRAY_TASK_ID"
export MATLAB_LOG_DIR="$SCRATCH/mathworks_cache/logs"
mkdir -p "$MATLAB_PREFDIR" "$MATLAB_LOG_DIR"

export CSC_STAGE=slice

cd "$CODE_DIR"
echo "Host: $(hostname)  |  Job: $SLURM_JOB_ID  |  Slice: $SLURM_ARRAY_TASK_ID  |  Started: $(date)"
matlab -nodisplay -nosplash -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); DogAnalysis"
echo "Finished: $(date)"
