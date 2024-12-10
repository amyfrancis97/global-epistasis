# Global Epistasis Models of Antibody-Antigen Interactions

## Introduction and Overview

This repository contains scripts associated with the Global Epistasis internship project. 

## Directory Structure
```bash
.
├── README.md                             # Project overview and usage instructions
├── config.sh                             # Configuration settings for the pipeline- **This must be updated for the models to work!**
├── main.sh                               # Main script used for executing pipeline
├── 1_download_structures.sh              # Downloads the Absolut! structure files for complexes in data/global_epistasis_cdrs_greater_11.txt
├── 2_get_11_mer.sh                       # Gets the best-binding 11-mer slide for the above CDR-H3 sequences, for sequences >11 amino acids.
├── 3_get_mutants.sh                      # Executes get_mutants.py to generate all single mutants from the above 11-mer, and a specified percentage of doubles and triples. Then this script obtains Absolut! binding affinities for all mutants and re-formats for MAVE-NN.
├── data/                          
│   └── global_epistasis...txt            # Contains Absolut! antigen complexes & CDR-H3 sequence for which the WT CDR-H3 is predicted to be >=11 amino acids long by AbYsis.
├── get_mutants_scripts/
│   ├── get_mutants.py                    # Generates all singles and a specified number of doubles and triple mutants.
│   └── reformat_absolut_epitope.py       # Generates two mutant datasets- one constrained to a single epitope, and another considering epitope switching.
├── latent_phenotype_wrapper.sh           # Wrapper script for executing latent_phenotype_job.sh for each antigen complex.
├── latent_phenotype_job.sh               # Executes sampling_scripts/get_latent_phenotype_info.py
├── sampling_wrapper.sh                   # Wrapper script for executing sampling_strategies.sh for each antigen complex.
├── sampling_strategies.sh                # Executes different sub-sampling scripts.
├── sampling_scripts/
│   ├── config.py                         # Contains directory paths for input data and results. 
│   ├── utils.py                          # Contains common MAVE-NN functions and models
│   ├── sample_data_all.py                # Script for randomly sub-sampling training data
│   └── sample_data_doubles_triples.py    # Script for randomly sub-sampling training data
├── visualisation/
│   └── plot_latent_phenotypes.py         # Script for plotting latent phenotypes vs observed phenotypes.
└── esm/
    ├── get_antibert_likelihoods.py       # Mutate sequences based on AntiBert likelihoods
    ├── get_esm_likelihoods.py            # Mutate sequences based on ESM likelihoods
    └──  get_iglm_likelihoods.py          # Mutate sequences based on IGLM likelihoods

```

## Setup

### **Anaconda Setup:**
Install dependencies using conda (update ```<envname>``` to the environment name of your choice e.g., 'CanDrivR-env'):

```bash
conda env create --name <envname> --file=global_epistasis.yml
conda activate <envname>
```
## Usage
1. **Update the config file to contain the paths in which this repository has been cloned.**
   
```bash
nano config.sh
```



