# 🧬 Global Epistasis Models for Antibody-Antigen Interactions 

## 📖 Introduction and Overview 

This repository contains scripts associated with the Global Epistasis internship project. 

## 🌳 Directory Structure
```bash
.
├── README.md                             # Project overview and usage instructions.
├── config.sh                             # Configuration settings for the pipeline- **This must be updated for the models to work!**
├── run_jobs.sh                           # Contains specifications for executing code on Slurm compute systems- **This must be updated to match the machine specifications!**
├── main.sh                               # Main script used for executing the pipeline.
├── 1_get_structures.sh                   # Downloads the Absolut! structure files for complexes in data/global_epistasis_cdrs_greater_11.txt.
├── 2_get_11_mer.sh                       # Gets the best-binding 11-mer slide for the above CDR-H3 sequences, for sequences >11 amino acids.
├── 3_get_mutants.sh                      # Executes get_mutants.py to generate all single mutants from the above 11-mer, and a specified percentage of doubles and triples. Then this script obtains Absolut! binding affinities for all mutants and re-formats for MAVE-NN.                   
├── 4_run_ablation.sh                     # Runs ablation analysis for all antigen complexes.
├── 5_get_latent_phenotype_plots.sh       # Trains a model for a specified number of doubles and triples for each antigen complex, saves data to plot latent phenotype and prediction plots.
├── ablation_analysis.sh                  # Executes Python script for ablation analysis.
├── data/                          
│   └── global_epistasis...txt            # Contains Absolut! antigen complexes & CDR-H3 sequence for which the WT CDR-H3 is predicted to be >=11 amino acids long by AbYsis.
├── get_mutants_scripts/
│   ├── get_mutants.py                    # Generates all singles and a specified number of doubles and triple mutants.
│   └── reformat_absolut_epitope.py       # Generates two mutant datasets- one constrained to a single epitope, and another considering epitope switching.
├── sampling_scripts/
│   ├── config.py                         # Contains directory paths for input data and results. 
│   ├── utils.py                          # Contains common MAVE-NN functions and models
│   ├── get_latent_phenotype_info.py      # Runs models for a specified number of double and triple mutants and saves data for plotting latent phenotype and prediction plots.
│   ├── sample_data_all.py                # Script for randomly sub-sampling training data
│   └── sample_data_doubles_triples.py    # Script for randomly sub-sampling numbers of doubles and triples.
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
### **Absolut! Setup:**

In order to run these scripts, Absolut! software must first be set up. For this analyis, we only require the ./AbsolutNoLib version of the software to be set up, which allows analysis on pre-computed structures. For more information on how to set up this software, please follow the documentation in the [Absolut! GitHub Repository](https://github.com/csi-greifflab/Absolut).


## 💻 Usage
1. **Update the config file to contain the paths in which this repository has been cloned.**
   
```bash
nano config.sh
```

2. **Update the job file to match compute specification.**
   
```bash
nano run_jobs.sh
```
Currently, job script uses slurm. Theis job scripts will need to be updated to match the specifications of the machine that is in use.

2. **Optional: Edit and execute the main.sh script.**
   
```bash
chmod +x main.sh
nano main.sh # Update to compute specifications.
./main.sh
```

The main script executes the following job scripts:

* 1_download_structures.sh
* 2_get_11_mer.sh
* 3_get_mutants.sh
* 4_run_ablation.sh
* 5_get_latent_phenotype_plots.sh

For the best-binding 11-mer CDR-H3 sequence (obtained in the previous script), this executes scripts to generate all possible single mutants, and a specified percentage of double and triple mutants for each antigen in 'data/global_epistasis_cdrs_greater_11.txt'. Then it executes another script to obtain predicted Asbolut! binding affinities for each of the mutant sequences and reformats for MAVE-NN models.

The default is for this script to obtain binding affinities for all singles, 50% of doubles, and 1% of triples, but these percentages can be specified as arguments within the main script.

## ℹ️ Data Details: **data/global_epistasis_cdrs_greater_11.txt**
This text file contains the names of the antigen structure files and the corresponding CDR-H3 sequence. To obtain these sequences, the heavy chain sequence was obtained from the antibody using NCBI, and then the CDR-H3 region was predicted using abYsis. Only antigen structure files where the corresponding antibody CDR-H3 region was predicted to be >11 amino acids long was retained (due to Absolut! predictions working on 11-mer sliding windows). 

## ℹ️ Script Details

### Step 1: 1_get_structures.sh ###
This script reads in the antigen structure file names from *data/global_epistasis_cdrs_greater_11.txt*. It then gets the filepath info using ./AbsolutNoLib and downloads the corresponding structure files. 

Execution:
```bash
sbatch run_jobs.sh 1_get_structures.sh
```

### Step 2: 2_get_11_mer.sh ###
Once the structure files have been downloaded from the Absolut! software (above), this script then takes the corresponding CDR-H3 sequence from *data/global_epistasis_cdrs_greater_11.txt* and runs them through the software to obtain the predicted binding affinities. *./AbsolutNoLib singleBinding* predicts the binding of each 11-mer slide for a given sequence. Therefore, for CDR-H3 sequences longer than 11 amino acids long, this script then selects the best-binding 11-mer slide for each antigen complex (e.g., the most negative).

Execution:
```bash
sbatch run_jobs.sh 2_get_11_mer.sh $pdb_code $sequence
```

Where $pdb_code is the antigen structure file and the $sequence is the CDR-H3 sequence of the antibody. These are extracted from *data/global_epistasis_cdrs_greater_11.txt*.

### Step 3: 3_get_mutants.sh ###
For the best-binding 11-mer CDR-H3 sequence (obtained in the previous script), this executes python scripts that generate all possible single mutants, and a specified percentage of double and triple mutants. The default is to generate all single mutants, 50% of possible double mutants, and 5% of triple mutants.

Once the mutant sequences have been randomly generated, the script then runs *./AbsolutNoLib repertoire* to obtain the binding affinities for each of the mutants, and finally generates a file that can be used as input for MAVE-NN software. For each antigen, the following output files are generated:

* trainval (all_mutants)     # Trainval file for all mutant sequences across all epitopes
* test (all_mutants)         # Test file for all mutant sequences across all epitopes
* trainval (single epitope)  # Trainval file filtered to only keep sequences that are predicted to bind to the same epitope as the wild-type sequence.
* test (single_epitope)      # Test file filtered to only keep sequences that are predicted to bind to the same epitope as the wild-type sequence.

Execution:
```bash
sbatch run_jobs.sh 3_get_mutants.sh --double 0.5 --triple 0.01 "$result_file"
```

Where *$result_file* is the resulting binding affinity file for each antigen that contains the best-binding 11-mer CDR-H3.

### Step 4: run_ablation.sh ###
This script runs MAVE-NN global epistasis models for different numbers of double and triples. The *ablation_analysis.sh* script must have two arguments passed to it:

```bash
./ablation_analysis.sh "all_mutants" "sample_all"
```

* The first argument can either be "all_mutants" or "single_epitope", depending on whether you want to allow epitope switching in Absolut! data.
* The second argument passed to the wrapper script is either "sample_all" or "sample_double_triples". "sample_all" randomly samples a specified total number of doubles and triples, whereas "sample_double_triples" allows you to specify exact numbers of doubles and exact numbers of triples to sample from the training dataset for modeling. Note that both options will always include all single mutants in the training data.

### Step 5: 5_get_latent_phenotype_plots.sh ###
This script runs a model for each antigen, for a specified number of double and triple mutant sequences. It then records and saves the data required for plotting latent phenotypes and prediction plots.

```bash
./5_get_latent_phenotype_plots.sh "all_mutants" 1000 1000 # Train model using all singles, 1000 doubles and 1000 triples
```


