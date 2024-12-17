#!/bin/bash

source config.sh 

# Set directory paths
SBATCH_SCRIPT="$SCRIPT_DIR/run_jobs.sh"

# Specify the filter (either "all_mutants" or "single_epitope")
FILTER=$1  # First argument passed to the wrapper script, e.g., "all_mutants" or "single_epitope"
SAMPLING_TYPE=$2 # Second argument passed to the wrapper script, e.g., "sample_all" or "sample_double_triples"

# Check if the filter is valid
if [ "$FILTER" != "all_mutants" ] && [ "$FILTER" != "single_epitope" ]; then
    echo "Invalid filter specified. Please use 'all_mutants' or 'single_epitope'."
    exit 1
fi

# Create an array of trainval files
trainval_files=(${DATA_DIR}/*${FILTER}*_trainval_df.txt)

# Check if trainval files exist
if [ ${#trainval_files[@]} -eq 0 ]; then
    echo "No trainval files found matching the filter '${FILTER}'."
    exit 1
fi

# Iterate over all trainval/test pairs
for trainval_file in "${trainval_files[@]}"; do
    # Extract the base name from the trainval file
    base_name=$(basename "$trainval_file" "_trainval_df.txt")

    # Look for the corresponding test file
    test_file="${DATA_DIR}/${base_name}_test_df.txt"
    
    if [ -f "$test_file" ]; then
        # Submit a job for this pair of files
        echo "Submitting job for $base_name"
        sbatch $SBATCH_SCRIPT ablation_analysis.sh "seed" $SAMPLING_TYPE "gpmap_type" "$trainval_file" "$test_file"
    else
        echo "No matching test file for $base_name"
    fi
done

