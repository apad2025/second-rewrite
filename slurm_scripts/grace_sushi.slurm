#!/bin/bash
#############################################################################
# grace_sushi.slurm           Process only the last SUSHI date
#
# Submit using the command:   sbatch grace_sushi.slurm
# Check progress:             squeue -u $USER   /   tail -f fwsep_first_*.out
#############################################################################
#SBATCH --job-name=fwsep_first
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --output=fwsep_first_%j.out
#SBATCH --error=fwsep_first_%j.err
### ---- EDIT THESE ---------------------------------------------------------
#SBATCH --account=ACCOUNT_XXXXXX
#SBATCH --partition=short
#############################################################################

set -euo pipefail

# ---- EDIT: MATLAB module --------------------------------------------------
# Run `module spider MATLAB` on Grace to see exact names, then set below.
MATLAB_MODULE="Matlab/R2024a"

# ---- Path to the folder that contains DogAnalysis.m -----------------------
CODE_DIR="/scratch/user/apad/GRMD/Images"

# ---- Environment ----------------------------------------------------------
module purge
module load "$MATLAB_MODULE"

# ---- Run ------------------------------------------------------------------
# DogAnalysis adds GRMD + hernando to the path itself, and loads DD.mat via
# DogInitialize, so DD.mat must already exist (run dogexplorer once first).
# To regenerate DD.mat in this job too, uncomment the next line:
# matlab -nodisplay -nosplash -batch "addpath(genpath('/scratch/user/apad/GRMD')); dogexplorer"

cd "$CODE_DIR"
echo "Host: $(hostname)  |  Job: $SLURM_JOB_ID  |  Started: $(date)"
matlab -nodisplay -nosplash -batch "DogAnalysis"
echo "Finished: $(date)"

