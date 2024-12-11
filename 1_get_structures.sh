#!/bin/bash

source config.sh

cd $ABSOLUT_DIR

# Iterate over each line of the file
while read -r line; do
    # Extract the first column (PDB code)
    pdb_code=$(echo "$line" | awk '{print $1}')
    
    # Run the AbsolutNoLib command with the pdb_code
    ./AbsolutNoLib info_filenames "$pdb_code" > "${pdb_code}_file_info.txt"
    
    # Use grep to find the curl command and execute it
    bash <(grep 'curl' "${pdb_code}_file_info.txt")
    
done < ${SCRIPT_DIR}/data/global_epistasis_cdrs_greater_11.txt

# After the loop, unzip all downloaded .zip files in the current directory
for zip_file in *.zip; do
    if [ -f "$zip_file" ]; then
        # Check if the file is a valid zip file
        if file "$zip_file" | grep -q "Zip archive data"; then
            # Unzip the file, automatically overwrite existing files
            unzip -o "$zip_file"
            # Remove the .zip file after unzipping
            rm "$zip_file"
        else
            echo "$zip_file is not a valid zip file."
        fi
    fi
done

