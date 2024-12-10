#!/bin/bash

# Set directory paths
DATA_DIR="/user/home/uw20204/global-epistasis/mavenn_pipeline/pipeline2/pipeline_scaling/data"
SBATCH_SCRIPT="/user/home/uw20204/global-epistasis/mavenn_pipeline/pipeline2/pipeline_scaling/latent_phenotype_job.sh"
OUTPUT_PATH="/user/home/uw20204/global-epistasis/mavenn_pipeline/pipeline2/pipeline_scaling/results"

# Arguments
FILTER=$1        # Filter type, e.g., "all_mutants" or "single_epitope"
NUM_DOUBLES=$2   # Number of doubles to sample
NUM_TRIPLES=$3   # Number of triples to sample

# Max number of jobs to submit at a time
MAX_JOBS=3

# Validate inputs
if [ "$FILTER" != "all_mutants" ] && [ "$FILTER" != "single_epitope" ]; then
    echo "Invalid filter. Use 'all_mutants' or 'single_epitope'."
    exit 1
fi

# Collect trainval files based on the filter
trainval_files=(${DATA_DIR}/*${FILTER}*_trainval_df.txt)

# Iterate over all trainval/test file pairs
for trainval_file in "${trainval_files[@]}"; do
    # Extract base name and locate the matching test file
    base_name=$(basename "$trainval_file" "_trainval_df.txt")
    test_file="${DATA_DIR}/${base_name}_test_df.txt"

    if [ -f "$test_file" ]; then
        echo "Submitting job for $base_name..."

        # Submit the job using sbatch
        sbatch "$SBATCH_SCRIPT" \
            "$base_name" \
            "sample_double_triples" \
            "additive" \
            "$trainval_file" \
            "$test_file" \
            "$NUM_DOUBLES" \
            "$NUM_TRIPLES" \
            "$OUTPUT_PATH"

        # Check the number of running jobs and wait if MAX_JOBS is reached
        while [ "$(squeue -u $USER | wc -l)" -ge $((MAX_JOBS + 1)) ]; do
            echo "Max number of jobs ($MAX_JOBS) reached. Waiting..."
            sleep 90  # Check again after 10 seconds
        done
    else
        echo "No matching test file for $base_name. Skipping."
    fi
done

