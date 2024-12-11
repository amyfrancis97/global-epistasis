#!/bin/bash

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
rm "${pdb_code}_results.tmp"
