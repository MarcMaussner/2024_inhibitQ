#!/bin/bash

# Set the current directory as the working directory
WORK_DIR=$(pwd)

# Function to check if a command was successful
check_success() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        exit 1
    fi
}

# Run CP2K in Docker
echo "Starting CP2K calculation..."
docker run --volume-driver local -v $WORK_DIR:/mnt --shm-size=1g --rm --user root cp2k/cp2k sh -c "umask 0000 && mpiexec -genv OMP_NUM_THREADS=1 -np 16 cp2k Al111_active_space.inp" > cp2k.log 2>&1 &

# Store the Docker process ID
DOCKER_PID=$!

# Wait for the socket file to be created
while [ ! -S $WORK_DIR/embedding_socket ]; do
    sleep 1
done

# Run the Python script
echo "Starting Python VQE calculation..."
python -u client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5 2>&1 | tee python_output.log

# Wait for the Docker process to finish
wait $DOCKER_PID
check_success "CP2K calculation"

# Ensure all files are readable and writable
chmod -R a+rw $WORK_DIR

echo "Calculations completed. Check cp2k.log and python_output.log for results."