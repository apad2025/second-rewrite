#!/bin/bash
#############################################################################
# grace_allslices_parfor.slurm    Separate all 50 slices in parallel (parfor)
#
# Submit using the command:   sbatch grace_allslices_parfor.slurm
# Check progress:             squeue -u $USER   /   tail -f fwsep_par_*.out
#
# Function_Bipolar_GC now uses as parfor
#############################################################################
#SBATCH --job-name=fwsep_par
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=180G
#SBATCH --output=fwsep_par_%j.out
#SBATCH --error=fwsep_par_%j.err
### ---- EDIT THESE ---------------------------------------------------------
#SBATCH --account=ACCOUNT_XXXXXX
#SBATCH --partition=medium
#############################################################################

set -euo pipefail

MATLAB_MODULE="Matlab/R2024a"
CODE_DIR="/scratch/user/apad/GRMD/Images"
mkdir -p /scratch/user/apad/residuals

module purge
module load "$MATLAB_MODULE"

export MATLAB_PREFDIR="$SCRATCH/mathworks_cache/prefdir"   # preferences dir
export MATLAB_LOG_DIR="$SCRATCH/mathworks_cache/logs"      # crash/diagnostic logs
mkdir -p "$MATLAB_PREFDIR" "$MATLAB_LOG_DIR"

cd "$CODE_DIR"
echo "Host: $(hostname)  |  Job: $SLURM_JOB_ID  |  Cores: $SLURM_CPUS_PER_TASK  |  Started: $(date)"
matlab -nodisplay -nosplash -batch "DogAnalysis"
echo "Finished: $(date)"
