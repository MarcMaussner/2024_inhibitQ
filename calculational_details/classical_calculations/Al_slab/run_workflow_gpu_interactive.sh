#!/bin/bash

# Set the CP2K tag and number of MPI processes
export CP2K_TAG="v2023.2"
export NUM_MPI_PROCESSES=16

# Function to run CP2K calculation on GPU
run_cp2k_gpu() {
    local input_file=$1
    echo "Running CP2K calculation for $input_file on GPU with $NUM_MPI_PROCESSES MPI processes"
    docker run -it --rm \
        --gpus all \
        --shm-size 32Gb \
        -v ${PWD}:/host_pwd \
        --workdir /host_pwd \
        nvcr.io/hpc/cp2k:${CP2K_TAG} \
        bash -c "export OMP_NUM_THREADS=1 && \
                 export CP2K_DATA_DIR=/usr/local/cp2k/data && \
                 export CUDA_VISIBLE_DEVICES=0 && \
                 mpirun -np $NUM_MPI_PROCESSES cp2k.psmp -i $input_file"
    sleep 5
}

# Run CP2K calculation on GPU
run_cp2k_gpu Al111.inp > cp2k_output.log 2>&1

echo "GPU Workflow completed. Check cp2k_output.log for details."