#!/bin/bash

#SBATCH --job-name=new_mutants_job          # Job name
#SBATCH --output=logs/new_mutants_%j.out    # Output and error log
#SBATCH --error=logs/new_mutants_%j.err     # Error log
##SBATCH --partition=gpu                     # Partition with GPU resources
##SBATCH --gres=gpu:1                        # Request 1 GPU
#SBATCH --cpus-per-task=4                   # Number of CPU cores per task
#SBATCH --mem=10G                           # Total memory limit
#SBATCH --time=24:00:00                     # Time limit: 24 hours
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --account=sscm013903
##SBATCH --chdir=/user/home/uw20204/global-epistasis/mavenn_pipeline/pipeline2/pipeline_scaling

source config.sh

# Function to print usage instructions
usage() {
    echo "Usage: $0 [--double <percentage_double>] [--triple <percentage_triple>] <input_file>"
    echo "Example: $0 --double 0.5 --triple 0.05 input_file.txt"
    exit 1
}

# Set default values for double and triple mutants
DEFAULT_PERC_DOUBLE_MUTANTS=0.05
DEFAULT_PERC_TRIPLE_MUTANTS=0.005

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --double) PERC_DOUBLE_MUTANTS="$2"; shift ;;
        --triple) PERC_TRIPLE_MUTANTS="$2"; shift ;;
        *) input_file="$1"; shift ;;
    esac
    shift
done

# Ensure input file is provided
if [[ -z "$input_file" ]]; then
    usage
fi

# Use default values if parameters are not provided
PERC_DOUBLE_MUTANTS="${PERC_DOUBLE_MUTANTS:-$DEFAULT_PERC_DOUBLE_MUTANTS}"
PERC_TRIPLE_MUTANTS="${PERC_TRIPLE_MUTANTS:-$DEFAULT_PERC_TRIPLE_MUTANTS}"

THREADS=50

# Process each line of the input file
while read -r ANTIGEN _ CDR_SEQ _ _ EPITOPE; do

    echo "Processing Antigen: $ANTIGEN"
    echo "CDR Sequence: $CDR_SEQ"
    echo "Epitope: $EPITOPE"

    # Navigate to the working directory
    cd $ABSOLUT_DIR || { echo "Error: Could not change to directory $ABSOLUT_DIR"; exit 1; }

    start_time=$(date +%s)

    # Step 1: Generate mutants using the Python script
    generate_mutants() {
        echo "Generating mutants with $PERC_DOUBLE_MUTANTS double mutants and $PERC_TRIPLE_MUTANTS triple mutants for $CDR_SEQ..."
        python "$SCRIPT_DIR/get_mutants.py" \
            --sequence "$CDR_SEQ" \
            --perc_double_mutants "$PERC_DOUBLE_MUTANTS" \
            --perc_triple_mutants "$PERC_TRIPLE_MUTANTS" \
            --filename "${ANTIGEN}_${PERC_DOUBLE_MUTANTS}_${PERC_TRIPLE_MUTANTS}" \
            --output_dir "$DATA_DIR"

        if [ $? -ne 0 ]; then
            echo "Error: Python script to generate mutants failed"
            exit 1
        else
            echo "Mutants generation completed successfully."
        fi
    }

    # Step 2: Run AbsolutNoLib with generated mutants
    run_absolut() {
        echo "Running AbsolutNoLib with generated mutants..."
        REPERTOIRE_FILE="${DATA_DIR}/${ANTIGEN}_${PERC_DOUBLE_MUTANTS}_${PERC_TRIPLE_MUTANTS}.txt"
        FINAL_BINDINGS_FILE="${DATA_DIR}/${PERC_DOUBLE_MUTANTS}_${PERC_TRIPLE_MUTANTS}_${ANTIGEN}FinalBindings_Process_1_Of_1.txt"

        ./AbsolutNoLib repertoire $ANTIGEN "$REPERTOIRE_FILE" $THREADS "${DATA_DIR}/${PERC_DOUBLE_MUTANTS}_${PERC_TRIPLE_MUTANTS}_"

        if [ $? -ne 0 ]; then
            echo "Error: AbsolutNoLib failed at mutant run"
            exit 1
        else
            echo "AbsolutNoLib execution completed successfully."
        fi

        # Validate that the final bindings file is two lines longer than the repertoire file
        echo "Validating the length of the final bindings file..."
        REPERTOIRE_LINE_COUNT=$(wc -l < "$REPERTOIRE_FILE")
        FINAL_BINDINGS_LINE_COUNT=$(wc -l < "$FINAL_BINDINGS_FILE")

        if [ $FINAL_BINDINGS_LINE_COUNT -ne $((REPERTOIRE_LINE_COUNT + 2)) ]; then
            echo "Warning: The final bindings file does not have the expected length."
            echo "Repertoire file lines: $REPERTOIRE_LINE_COUNT"
            echo "Final bindings file lines: $FINAL_BINDINGS_LINE_COUNT"
        else
            echo "Validation successful: Final bindings file is two lines longer than the repertoire file."
        fi
    }

    # Step 3: Generate MAVENN input from AbsolutNoLib output
    generate_mavenn_input() {
        echo "Generating MAVENN input from AbsolutNoLib output..."
        python -u "$SCRIPT_DIR/reformat_absolut_epitope.py" \
            --absolut_final_bindings_path "${DATA_DIR}/${PERC_DOUBLE_MUTANTS}_${PERC_TRIPLE_MUTANTS}_${ANTIGEN}FinalBindings_Process_1_Of_1.txt" \
            --mutant_file_path_with_dist "${DATA_DIR}/${ANTIGEN}_${PERC_DOUBLE_MUTANTS}_${PERC_TRIPLE_MUTANTS}_with_dist.txt" \
            --output_dir "$DATA_DIR" \
            --cdr_selection "$CDR_SEQ" \
            --antigen "$ANTIGEN" \
            --percentage_prefix "${PERC_DOUBLE_MUTANTS}_${PERC_TRIPLE_MUTANTS}" \
            --filter_epitope "$EPITOPE"

        if [ $? -ne 0 ]; then
            echo "Error: MAVENN input generation failed"
            exit 1
        else
            echo "MAVENN input generation completed successfully."
        fi
    }

    # Main execution steps
    generate_mutants
    run_absolut
    generate_mavenn_input

    # Remove all files except those ending with "trainval_df.txt" or "test_df.txt"
    echo "Cleaning up intermediate files..."
    cd "$DATA_DIR"
    rm *${PERC_DOUBLE_MUTANTS}_${PERC_TRIPLE_MUTANTS}_Temp*
    echo "Cleanup completed."

    finish_time=$(date +%s)
    echo "Total time taken for $ANTIGEN: $((finish_time - start_time)) seconds."

done < "$input_file"


