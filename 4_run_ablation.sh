#!/bin/bash

source config.sh 

# Set directory paths
SBATCH_SCRIPT="$SCRIPT_DIR/run_jobs.sh"

# Specify the filter (either "all_mutants" or "single_epitope")
FILTER=$1  # First argument passed to the wrapper script, e.g., "all_mutants" or "single_epitope"
SAMPLING_TYPE=$2 # Second argument passed to the wrapper script, e.g., "sample_all" or "sample_double_triples"
BATCH=$3  # Third argument passed to the wrapper script, e.g., "first_half" or "second_half"

# Check if the filter is valid
if [ "$FILTER" != "all_mutants" ] && [ "$FILTER" != "single_epitope" ]; then
    echo "Invalid filter specified. Please use 'all_mutants' or 'single_epitope'."
    exit 1
fi

# Create an array of trainval files
trainval_files=(${DATA_DIR}/*${FILTER}*_trainval_df.txt)

# Count total number of files
total_files=${#trainval_files[@]}
half=$((total_files / 2))

# Determine the range of indices based on the batch
if [ "$BATCH" == "first_half" ]; then
    start_idx=0
    end_idx=$half
elif [ "$BATCH" == "second_half" ]; then
    start_idx=$half
    end_idx=$total_files
else
    echo "Invalid batch specified. Please use 'first_half' or 'second_half'."
    exit 1
fi

# Iterate over the selected range of trainval/test pairs
for ((i=start_idx; i<end_idx; i++)); do
    trainval_file=${trainval_files[$i]}
    
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

