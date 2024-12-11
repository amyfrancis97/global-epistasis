#!/bin/bash

#SBATCH --job-name=global_epistasis            # Job name
#SBATCH --output=logs/jobs_%j.out              # Output and error log
#SBATCH --error=logs/jobs_%j.err               # Error log
#SBATCH --partition=gpu                        # Partition with GPU resources
#SBATCH --gres=gpu:1                           # Request 1 GPU
#SBATCH --cpus-per-task=4                      # Number of CPU cores per task
#SBATCH --mem=2G                              # Total memory limit
#SBATCH --time=5:00:00                        # Time limit: 24 hours
#SBATCH --ntasks=1                             # Run a single task
#SBATCH --account=sscm013903
#SBATCH --chdir=/user/home/uw20204/global-epistasis/mavenn_pipeline/pipeline2/pipeline_scaling

# Load configuration if needed
source config.sh

# Get the target script and shift arguments
TARGET_SCRIPT=$1
shift

# Run the target script with the remaining arguments
bash "$TARGET_SCRIPT" "$@"
