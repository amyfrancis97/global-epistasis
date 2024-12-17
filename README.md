# Global Epistasis Models for Antibody-Antigen Interactions 

## Introduction and Overview 

This repository contains scripts for the Global Epistasis Project, which focuses on predicting antibody-antigen interactions for double and triple mutant antibody sequences. The pipeline generates mutant datasets for a given CDR-H3 region, computes their binding affinities using Absolut!, and trains global epistasis models (MAVE-NN) to predict binding affinities of unseen double and triple mutants.

##  Directory Structure
```bash
.
├── README.md                                # Project instructions
├── config.sh                                # Configuration file for the pipeline (must be updated!)
├── run_jobs.sh                              # Job specifications for Slurm compute systems
├── main.sh                                  # Main script to execute the entire pipeline
├── 1_get_structures.sh                      # Downloads Absolut! structure files for specified complexes
├── 2_get_11_mer.sh                          # Finds best-binding 11-mer slides for CDR-H3 sequences
├── 3_get_mutants.sh                         # Generates mutants and computes binding affinities
├── 4_run_ablation.sh                        # Runs ablation analysis for antigen complexes
├── 5_get_latent_phenotype_plots.sh          # Trains models and generates phenotype plots
├── ablation_analysis.sh                     # Executes ablation analysis Python script
├── data/                                    # Input files and data
│   └── global_epistasis_cdrs_greater_11.txt # List of antigens and CDR-H3 sequences >11 amino acids
├── get_mutants_scripts/
│   ├── get_mutants.py                       # Generates mutant datasets (singles, doubles, triples)
│   └── reformat_absolut.py                  # Reformats data for MAVE-NN models.
├── sampling_scripts/
│   ├── config.py                            # Directory paths for input and output
│   ├── utils.py                             # Shared functions for MAVE-NN models
│   ├── get_latent_phenotype_info.py         # Saves latent phenotype data for plotting
│   ├── sample_data_all.py                   # Sub-samples training data
│   └── sample_data_doubles_triples.py       # Sub-samples double and triple mutants
├── visualisation/
│   ├── config.py                            # Config file for visualisation scripts and notebooks
│   ├── ablation_analysis.ipynb              # Notebook for visualising ablation analysis
│   ├── compare_epitope_switching.ipynb      # Notebook for comparing results with and without epitope switching
│   ├── mavenn_visualisation_modules.py      # Contains all visualisation functions
│   ├── plot_affinity_distributions.ipynb    # Visualise binding affinity distributions for mutants from a single antigen complex
│   └── plot_latent_phenotypes.py            # Plots observed vs predicted phenotypes
└── antibody_language_models/
    ├── get_antibert_likelihoods.py          # Generates mutations using AntiBERT likelihoods
    ├── get_esm_likelihoods.py               # Generates mutations using ESM likelihoods
    └── get_iglm_likelihoods.py              # Generates mutations using IGLM likelihoods

```

##  Setup

### **Anaconda Setup:**
Install dependencies using conda:

```bash
conda env create --name <envname> --file=global_epistasis.yml
conda activate <envname>
```

Replace `<envname>` with your desired environment name, e.g., `global_epistasis`.

### **Absolut! Setup:**
To run these scripts, you need to set up Absolut! software. For this project, only the `./AbsolutNoLib` version is required, which works with pre-computed structures. Follow the documentation in the  [Absolut! GitHub Repository](https://github.com/csi-greifflab/Absolut).


## Usage
1. **Update the configuration file** with repository paths:
   
```bash
nano config.sh
```

2. **Update the job script** to match your compute specifications (e.g., Slurm settings):
   
```bash
nano run_jobs.sh
```

2. **Make the main script executable and customise if necessary**:
   
```bash
chmod +x main.sh
nano main.sh # Update to compute specifications.
./main.sh
```
The main.sh script executes the following jobs:

* `1_get_structures.sh`
* `2_get_11_mer.sh`
* `3_get_mutants.sh`
* `4_run_ablation.sh`
* `5_get_latent_phenotype_plots.sh`

## **Data Details**: 
`data/global_epistasis_cdrs_greater_11.txt`
This file contains antigen structure file names and corresponding antibody CDR-H3 sequences predicted to be ≥11 amino acids long. The sequences were:

Extracted from antibody heavy chains using NCBI.
1. Analysed for CDR-H3 regions using abYsis.
2. Antigens with CDR-H3 lengths >11 amino acids were retained, as Absolut! predictions require 11-mer sliding windows.

## Script Details

#### Step 1: `1_get_structures.sh`
Downloads antigen structure files listed in `data/global_epistasis_cdrs_greater_11.txt` using **Absolut!**.

**Execution**:
```bash
sbatch run_jobs.sh 1_get_structures.sh
```

#### Step 2: `2_get_11_mer.sh`
Finds the best-binding 11-mer slide for each CDR-H3 sequence.

Execution:
```bash
sbatch run_jobs.sh 2_get_11_mer.sh <pdb_code> <sequence>
```
Replace `<pdb_code>` with the antigen file name and `<sequence>` with the CDR-H3 sequence.

#### Step 3: `3_get_mutants.sh`
Generates mutants and computes **Absolut!** binding affinities:

* Singles: All possible single mutants
* Doubles/Triples: Specified percentages (default: 50% doubles, 5% triples)

**Execution:**
```bash
sbatch run_jobs.sh 3_get_mutants.sh --double 0.5 --triple 0.01 <result_file>
```

#### Step 4:`4_run_ablation.sh`
Runs MAVE-NN global epistasis models.

```bash
./ablation_analysis.sh "all_mutants" "sample_all"
```

* `"all_mutants"`: Includes all mutants. Use `"single_epitope"` for epitope-constrained analysis.
* `"sample_all"`: Randomly samples doubles/triples. Use `"sample_double_triples"` for specific sampling.

#### Step 5: `5_get_latent_phenotype_plots.sh`
Trains models for specified numbers of doubles/triples and plots latent phenotypes.

```bash
./5_get_latent_phenotype_plots.sh "all_mutants" 1000 1000
```
Here, 1000 doubles and 1000 triples are sampled for training.

### Notes
* Ensure Absolut! software is set up correctly before executing the scripts.
* Update paths and Slurm configurations in config.sh and run_jobs.sh.
* Percentages for double and triple mutants can be adjusted in main.sh.
