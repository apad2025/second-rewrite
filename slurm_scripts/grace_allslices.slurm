#!/bin/bash
#############################################################################
# grace_allslices.slurm       Separate all 50 slices of last SUSHI date
#
# Submit using the command:   sbatch grace_allslices.slurm
# Check progress:             squeue -u $USER   /   tail -f fwsep_all_*.out
#
# Timing note: slice 12 took ~712 s. The slice loop is serial, so
# 50 slices ~= 50 * 712 s ~= 9.9 h. medium allows up to 1 day (24 h),
# so 14 h requested here leaves ample margin.
#############################################################################
#SBATCH --job-name=fwsep_all
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --output=fwsep_all_%j.out
#SBATCH --error=fwsep_all_%j.err
### ---- EDIT THESE ---------------------------------------------------------
#SBATCH --account=ACCOUNT_XXXXXX
#SBATCH --partition=medium
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
