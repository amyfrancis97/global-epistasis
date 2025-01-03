#!/bin/bash
start_time=$(date +%s)

# Receive arguments from sbatch
SEQ_SOURCE=$1        # First argument passed, e.g. "random" or "seed"
SCRIPT_TO_RUN=$2     # Second argument passed, either "active_learning", "sample_all", or "sample_double_triples"
GPMAP_TYPE=$3
TRAINVAL_DF=$4       # Trainval file path passed from wrapper script
TEST_DF=$5           # Test file path passed from wrapper script

source config.sh

# Activate environment and set paths
cd $WORK_DIR

cd $SCRIPT_DIR/sampling_scripts

# Define data paths

# Check which script to run based on the argument passed
if [ "$SCRIPT_TO_RUN" == "sample_all" ]; then
    echo "Running Sample All Script"
    python -u sample_data_all.py \
        --name $SEQ_SOURCE \
        --trainval_file ${TRAINVAL_DF} \
        --test_file ${TEST_DF} \
        --output_path ${OUTPUT_PATH}/${SEQ_SOURCE} \
        --num_samples_list "2-48-2,50-3000-25,4000-17000-2000" \
        --gpmap_types "additive"
##        --gpmap_types "additive,neighbor,pairwise,blackbox"

elif [ "$SCRIPT_TO_RUN" == "sample_double_triples" ]; then
    echo "Running Sample Double and Triple Mutants"
    python -u sample_data_doubles_triples.py \
        --name $SEQ_SOURCE \
        --trainval_file ${TRAINVAL_DF} \
        --test_file ${TEST_DF} \
        --output_path ${OUTPUT_PATH}/${SEQ_SOURCE} \
        --gpmap_types "additive" \
##        --gpmap_types "additive,neighbor,pairwise,blackbox" \
        --range_to_sample "1,2-48-2,50-2000-50,3000-17000-2000"

else
    echo "Invalid script name provided. Please use 'active_learning', 'sample_all', or 'sample_double_triples'."
    exit 1
fi

finish_time=$(date +%s)
echo "Total time taken for $SEQ_SOURCE: $((finish_time - start_time)) seconds."

