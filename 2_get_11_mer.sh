#!/bin/bash
#SBATCH --job-name=get_structures_job          # Job name
#SBATCH --output=logs/get_structures_%j.out         # Output and error log
#SBATCH --error=logs/get_structures_%j.err          # Error log
#SBATCH --partition=gpu                    # Partition with GPU resources
#SBATCH --gres=gpu:1                       # Request 1 GPU
#SBATCH --cpus-per-task=4                   # Number of CPU cores per task
#SBATCH --mem=10G                           # Total memory limit
#SBATCH --time=24:00:00                     # Time limit: 24 hours
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --account=sscm013903
#SBATCH --chdir=/user/home/uw20204/global-epistasis/mavenn_pipeline/pipeline2/pipeline_scaling

source config.sh

pdb_code=$1
sequence=$2

echo $pdb_code
echo $sequence

cd $ABSOLUT_DIR

./AbsolutNoLib singleBinding "$pdb_code" "$sequence" > "${pdb_code}_results.tmp"

wait

# Remove the headers and keep the structures
awk '/=== Structures ready! ===/ {found=1; next} found'  "${pdb_code}_results.tmp" | sed '1d' > "${pdb_code}_results.txt"
awk -v pdb_code="$pdb_code" '{print pdb_code "\t" $0}' "${pdb_code}_results.txt" > "${pdb_code}_results.tmp"
awk -F"\t" 'NR == 1 {min=$3; line=$0} $3 < min {min=$3; line=$0} END {print line}' "${pdb_code}_results.tmp" > "${pdb_code}_results.txt"

