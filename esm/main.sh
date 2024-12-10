#!/bin/bash

#SBATCH --job-name=esm_jobs     # Job name
#SBATCH --output=logs/esm%j.out    # Output and error log
#SBATCH --error=logs/esm%j.err     # Error log
##SBATCH --partition=gpu                 # Partition with GPU resources
##SBATCH --gres=gpu:1                    # Request 1 GPU
#SBATCH --mem=12G                       # Total memory limit
#SBATCH --time=24:00:00                 # Time limit: 48 hours
#SBATCH --account=sscm013903
#SBATCH --chdir=/user/home/uw20204/global-epistasis/mavenn_pipeline/pipeline2/pipeline_scaling/esm

source /user/home/uw20204/.bashrc

conda activate global_epistasis

#python -u get_esm_likelihoods.py

python -u get_antibert_likelihoods.py

#python -u get_iglm_likelihoods.py
