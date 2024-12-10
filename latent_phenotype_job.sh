#!/bin/bash

#SBATCH --job-name=sample_antigen     # Job name
#SBATCH --output=logs/%x_%j.out       # Output log
#SBATCH --error=logs/%x_%j.err        # Error log
#SBATCH --partition=gpu               # GPU partition
#SBATCH --gres=gpu:1                  # Request 1 GPU
#SBATCH --mem=64G                     # Memory allocation
#SBATCH --time=48:00:00               # Max time for job
#SBATCH --account=sscm013903          # Account name
#SBATCH --chdir=/user/home/uw20204/global-epistasis/mavenn_pipeline/pipeline2/pipeline_scaling

# Load environment
source /user/home/uw20204/.bashrc
conda activate global_epistasis

# Receive arguments from submission script
BASE_NAME=$1
SCRIPT_TYPE=$2
GPMAP_TYPES=$3
TRAINVAL_FILE=$4
TEST_FILE=$5
NUM_DOUBLES=$6
NUM_TRIPLES=$7
OUTPUT_PATH=$8

# Navigate to script directory
cd /user/home/uw20204/global-epistasis/mavenn_pipeline/pipeline2/pipeline_scaling/sampling_scripts

# Run the sampling Python script
python get_latent_phenotype_info.py \
    --name "$BASE_NAME" \
    --trainval_file "$TRAINVAL_FILE" \
    --test_file "$TEST_FILE" \
    --output_path "$OUTPUT_PATH" \
    --gpmap_types "$GPMAP_TYPES" \
    --num_doubles "$NUM_DOUBLES" \
    --num_triples "$NUM_TRIPLES"

