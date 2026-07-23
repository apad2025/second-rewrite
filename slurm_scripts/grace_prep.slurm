#!/bin/bash
#############################################################################
# grace_prep.slurm           Stage 1/3 of the slice-wise separation pipeline
#
# Preprocesses the dataset (zero-fill, mask, trim), saves RAWzip.mat, and
# caches the volume-wide SNR threshold that every slice job has to share.
#
# Submit using the command:   sbatch grace_prep.slurm
# Then, in order:             sbatch grace_slices.slurm
#                             sbatch grace_combine.slurm
#
# Skip straight to grace_slices.slurm if RAWzip.mat and SNRThreshold.mat
# already exist for this dog & date.
#############################################################################
#SBATCH --job-name=fwsep_prep
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --output=fwsep_prep_%j.out
#SBATCH --error=fwsep_prep_%j.err
### ---- EDIT THESE ---------------------------------------------------------
#SBATCH --account=ACCOUNT_XXXXXX
#SBATCH --partition=medium
#############################################################################

set -euo pipefail

MATLAB_MODULE="Matlab/R2024a"
CODE_DIR="/scratch/user/apad/GRMD/Images"

module purge
module load "$MATLAB_MODULE"

export MATLAB_PREFDIR="$SCRATCH/mathworks_cache/prefdir"   # preferences dir
export MATLAB_LOG_DIR="$SCRATCH/mathworks_cache/logs"      # crash/diagnostic logs
mkdir -p "$MATLAB_PREFDIR" "$MATLAB_LOG_DIR"

export CSC_STAGE=prep

cd "$CODE_DIR"
echo "Host: $(hostname)  |  Job: $SLURM_JOB_ID  |  Stage: $CSC_STAGE  |  Started: $(date)"
matlab -nodisplay -nosplash -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); DogAnalysis"
echo "Finished: $(date)"
