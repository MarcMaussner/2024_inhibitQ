#!/bin/bash

# Function to run CP2K calculation
run_cp2k() {
    local input_file=$1
    local cp2k_version="2024.1"  # Updated to the version you have
    echo "Running CP2K calculation for $input_file"
    docker run --platform linux/amd64 -v $PWD:/mnt --shm-size=1g -it --rm --name cp2k_${cp2k_version} cp2k/cp2k:${cp2k_version} mpirun -np 1 cp2k -i /mnt/$input_file
    sleep 5
}

# Run CP2K calculation
run_cp2k MgO.inp

# Run analysis
python3 analysis.py

echo "Workflow completed."