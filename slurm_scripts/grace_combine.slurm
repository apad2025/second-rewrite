#!/bin/bash
#############################################################################
# grace_combine.slurm        Stage 3/3 of the slice-wise separation pipeline
#
# Stacks the per-slice saves into the full corrected dataset and writes it out
# under the usual FW name, so that it is interchangeable with the file the
# whole-volume run produces. Errors out listing any slices that are missing.
#
# Submit using the command:   sbatch grace_combine.slurm
#
# Background field removal & dipole inversion are 3D and are not run here;
# run DogAnalysis without CSC_STAGE once this file exists.
#############################################################################
#SBATCH --job-name=fwsep_comb
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --output=fwsep_comb_%j.out
#SBATCH --error=fwsep_comb_%j.err
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

export CSC_STAGE=combine

cd "$CODE_DIR"
echo "Host: $(hostname)  |  Job: $SLURM_JOB_ID  |  Stage: $CSC_STAGE  |  Started: $(date)"
matlab -nodisplay -nosplash -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); DogAnalysis"
echo "Finished: $(date)"
