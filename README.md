# 🧬 Global Epistasis Models of Antibody-Antigen Interactions 

## 📖 Introduction and Overview 

This repository contains scripts associated with the Global Epistasis internship project. 

## 🌳 Directory Structure
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

## 🛠️ Setup

### 🐍 **Anaconda Setup:**
Install dependencies using conda (update ```<envname>``` to the environment name of your choice e.g., 'global_epistasis'):

```bash
conda env create --name <envname> --file=global_epistasis.yml
conda activate <envname>
```
## 💻 Usage
1. **Update the config file to contain the paths in which this repository has been cloned.**
   
```bash
nano config.sh
```

2. **Edit and execute the main.sh script.**
   
```bash
./main.sh
```
The main script executes the following job scripts:

```bash
1_download_structures.sh
```
This script downloads the Absolut! structure files for the antigen complexes listed in 'data/global_epistasis_cdrs_greater_11.txt'

```bash
2_get_11_mer.sh
```
This script then obtains the binding affinities for each of the structures from Absolut!. If the CDR-H3 sequence is greater than 11 amino acids, this script obtains the binding affinity for each 11-mer slide, and returns the 11-mer CDR sequence that results in the most negative binding affinity.

```bash
3_get_mutants.sh
```

For the best-binding 11-mer CDR-H3 sequence (obtained in the previous script), this executes scripts to generate all possible single mutants, and a specified percentage of double and triple mutants. Then it executes another script to obtain predicted Asbolut! binding affinities for each of the mutant sequences and reformats for MAVE-NN models.

The default is for this script to obtain binding affinities for all singles, 50% of doubles, and 1% of triples, but these percentages can be specified as arguments within the main script:

```bash
sbatch 3_get_mutants.sh --double 0.5 --triple 0.01 "$result_file"
```

Then the script exectutes:

```bash
sampling_wrapper.sh
```

This job script runs MAVE-NN global epistasis models for different numbers of double and triples. sampling_wrapper.sh must have 3 arguments passed to it:

```bash
./sampling_wrapper.sh "all_mutants" "sample_all" "first_half"
```

The first argument passed to the wrapper script can either be "all_mutants" or "single_epitope", depending on whether you want to allow epitope switching in Absolut! data.
The second argument passed to the wrapper script is either "sample_all" or "sample_double_triples". "sample_all" randomly samples a specified total number of doubles and triples, whereas "sample_double_triples" allows you to specify exact numbers of doubles and exact numbers of triples to sample from the training dataset for modeling.
The third argument passed to the wrapper script breaks down the antigen complexes into two batches, implemented due to computational constraints e.g., "first_half" or "second_half".
