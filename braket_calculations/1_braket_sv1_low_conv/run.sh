#!/bin/bash

# Set AWS credentials and region
export AWS_DEFAULT_REGION=us-east-1  #eu-north-1
export AWS_ACCESS_KEY_ID=AKIAQLVQQSYW27WSICDH
export AWS_SECRET_ACCESS_KEY=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

# Set the current directory as the working directory
WORK_DIR=$(pwd)

# Function to check if a command was successful
check_success() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed" >> $WORK_DIR/run.log
        exit 1
    fi
}

# Run CP2K in Docker
echo "Starting CP2K calculation..." >> $WORK_DIR/run.log
docker run --volume-driver local -v $WORK_DIR:/mnt --shm-size=1g --rm --user root cp2k/cp2k sh -c "umask 0000 && mpiexec -genv OMP_NUM_THREADS=2 -np 104 cp2k Al111_active_space.inp" > $WORK_DIR/cp2k.log 2>&1 &

# Store the Docker process ID
DOCKER_PID=$!

# Wait for the socket file to be created
while [ ! -S $WORK_DIR/embedding_socket ]; do
   sleep 1
done

# Run the Python script
echo "Starting Python VQE calculation..." >> $WORK_DIR/run.log
python -u client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5 --braket > $WORK_DIR/python_output.log 2>&1 &

#PYTHON_PID=$!

# Wait for both processes to finish
wait $DOCKER_PID
check_success "CP2K calculation"

#wait $PYTHON_PID
#check_success "Python VQE calculation"

# Ensure all files are readable and writable
chmod -R a+rw $WORK_DIR

echo "Calculations completed. Check cp2k.log and python_output.log for results." >> $WORK_DIR/run.log
