#!/bin/bash

# Initialize counter
count=0
max_files=25

# Loop through each zip file and extract the specified file
for zip_file in /data/leuven/358/vsc35887/master_thesis/Master-thesis-Wout-Vanderheyden/AF3_data/H1N1/*.zip; do
    # Break loop if 15 files have been processed
    if [[ $count -ge $max_files ]]; then
        break
    fi
    
    unzip -j "$zip_file" "*model_0*" -d /data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/1_cif_files
    
    # Increment counter
    ((count++))
done


# to run: chmod +x unzip.sh
# to run: ./unzip.sh

