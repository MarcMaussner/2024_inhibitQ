#!/bin/bash

# Set the CP2K tag and number of MPI processes
export CP2K_TAG="v2023.2"
export NUM_MPI_PROCESSES=16  # Adjust this based on your system's capabilities

# Function to run CP2K calculation on GPU
run_cp2k_gpu() {
    local input_file=$1
    echo "Running CP2K calculation for $input_file on GPU with $NUM_MPI_PROCESSES MPI processes"
    docker run -it --privileged --rm \
        --gpus all \
        --shm-size 32Gb \
        -v ${PWD}:/host_pwd \
        --workdir /host_pwd \
        nvcr.io/hpc/cp2k:${CP2K_TAG} \
        mpirun --bind-to none -n ${NUM_MPI_PROCESSES} cp2k.psmp -i $input_file
    sleep 5
}

# Run CP2K calculation on GPU
run_cp2k_gpu MgO.inp

# Run analysis
python3 analysis.py

echo "GPU Workflow completed."

#./run_workflow_gpu.sh | tee -a run_workflow.out
#time ./run_workflow_gpu.sh >> run_workflow.out 2>&1