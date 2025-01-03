#!/bin/bash

# Set the directory path for Absolut
source config.sh 

# Download the Absolut! structure files
sbatch run_jobs.sh 1_get_structures.sh

cd $ABSOLUT_DIR

# Obtain the best-binding 11-mer slide of the CDR sequence for each antibody
# Obtains the corresponding epitope for the selected CDR sequence
cat ${SCRIPT_DIR}/data/global_epistasis_cdrs_greater_11.txt | while read -r line; do
    # Extract the first column (PDB code) and second column (amino acid sequence)
    pdb_code=$(echo "$line" | awk '{print $1}')
    sequence=$(echo "$line" | awk '{print $2}')
    echo $pdb_code
    echo $sequence
    cd $SCRIPT_DIR
    #  Create a job script to submit
    sbatch ${SCRIPT_DIR}/run_jobs.sh 2_get_11_mer.sh $pdb_code $sequence;
done

# For each of the best-binding 11-mers selected for each antibody.
# Generate all single mutants, and a specified proportion of double and triple mutants.
# Outputs train_val and test data into the /data directory.

# Loop through each file that matches the pattern *results.txt
# Gets mutants for each antigen-antibody complex
for result_file in "$ABSOLUT_DIR"/*results.txt; do
    # Submit the job using sbatch for each result file
    echo "Submitting job for $result_file"
    sbatch ${SCRIPT_DIR}/run_jobs.sh 3_get_mutants.sh 0.00005 0.000001 "$result_file"
done

# Run models for different samples of data including epitope switching
./4_run_ablation.sh "all_mutants" "sample_all"
./4_run_ablation.sh "all_mutants" "sample_double_triples"

# Repeat but restrict to only mutations predicted to bind to the same epitope as the wild type CDR sequence
./4_run_ablation.sh "single_epitope" "sample_all" 
./4_run_ablation.sh "single_epitope" "sample_double_triples"

# Get latent phenotype mappings for a given number of doubles and triples
# This script saves all of the data required for generating latent phenotype plots into "results/plot_data/"
./5_get_latent_phenotype_plots.sh "all_mutants" 1000 1000 # Train model using all singles, 1000 doubles and 1000 triples
./5_get_latent_phenotype_plots.sh "single_epitope" 1000 1000 # Train model using all singles, 1000 doubles and 1000 triples

