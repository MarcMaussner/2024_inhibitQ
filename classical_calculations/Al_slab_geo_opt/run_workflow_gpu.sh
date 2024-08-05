#!/bin/bash

# Set the CP2K tag and number of MPI processes
export CP2K_TAG="v2023.2"
export NUM_MPI_PROCESSES=4  # You can experiment with this number

# Function to run CP2K calculation on GPU
run_cp2k_gpu() {
    local input_file=$1
    local output_file=$2
    echo "Running CP2K calculation for $input_file on GPU with $NUM_MPI_PROCESSES MPI processes"
    docker run -d --rm \
        --gpus all \
        --shm-size 32Gb \
        -v ${PWD}:/host_pwd \
        --workdir /host_pwd \
        nvcr.io/hpc/cp2k:${CP2K_TAG} \
        bash -c "export OMP_NUM_THREADS=1 && \
                 export CP2K_DATA_DIR=/usr/local/cp2k/data && \
                 export CUDA_VISIBLE_DEVICES=0 && \
                 mpirun -np $NUM_MPI_PROCESSES cp2k.psmp -i $input_file > $output_file 2>&1"
}

# Run CP2K calculation on GPU in the background
INPUT_FILE="Al111.inp"
OUTPUT_FILE="cp2k_output.log"
run_cp2k_gpu $INPUT_FILE $OUTPUT_FILE

# Get the Docker container ID
CONTAINER_ID=$(docker ps -q --filter ancestor=nvcr.io/hpc/cp2k:${CP2K_TAG})

echo "CP2K calculation started in background."
echo "To check the progress, use: docker logs -f $CONTAINER_ID"
echo "To stop the calculation, use: docker stop $CONTAINER_ID"
echo "Output will be saved to: $OUTPUT_FILE"

#docker logs -f <CONTAINER_ID>
#docker stop <CONTAINER_ID>
#docker ps
#tail -f cp2k_output.log