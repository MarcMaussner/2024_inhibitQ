#!/bin/bash

# Function to run CP2K calculation on M1 Mac (kept for compatibility)
run_cp2k_m1() {
    local input_file=$1
    local cp2k_version="2024.1"
    echo "Running CP2K calculation for $input_file on M1 Mac"
    docker run --platform linux/amd64 -v $PWD:/mnt --shm-size=1g -it --rm \
        --name cp2k_${cp2k_version} cp2k/cp2k:${cp2k_version} \
        mpirun -np 1 cp2k -i /mnt/$input_file
    sleep 5
}

# Function to run CP2K calculation on Linux on CPU only
run_cp2k_linux() {
    local input_file=$1
    local cp2k_version="2024.1_openmpi_generic_psmp"
    local num_physical_cores=18  # Your CPU has 18 physical cores
    local num_mpi_processes=$num_physical_cores
    local num_omp_threads=2  # Use 2 OpenMP threads per MPI process (utilizing hyperthreading)

    echo "Running CP2K calculation for $input_file on Linux"
    docker run -v $PWD:/mnt --shm-size=32g -it --rm \
        --name cp2k_${cp2k_version} cp2k/cp2k:${cp2k_version} \
        bash -c "export OMP_NUM_THREADS=${num_omp_threads} && \
                 mpirun -np ${num_mpi_processes} --map-by socket:PE=${num_omp_threads} \
                 cp2k.psmp -i /mnt/$input_file"
    sleep 5
}

# Detect the operating system
if [[ "$(uname)" == "Darwin" ]]; then
    # macOS (M1)
    run_cp2k_m1 MgO.inp
else
    # Linux
    run_cp2k_linux MgO.inp
fi

# Run analysis
python3 analysis.py

echo "Workflow completed."