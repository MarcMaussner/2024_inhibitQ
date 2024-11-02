#!/bin/bash

# Set AWS credentials and region
export AWS_DEFAULT_REGION=us-east-1  #eu-north-1

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
docker run --volume-driver local -v $WORK_DIR:/mnt --shm-size=1g --rm --user root cp2k/cp2k sh -c "umask 0000 && mpiexec -genv OMP_NUM_THREADS=1 -np 48 cp2k Al111_active_space.inp" > $WORK_DIR/cp2k.log 2>&1

# Remove the & to run in foreground
check_success "CP2K calculation"

# Ensure all files are readable and writable
chmod -R a+rw $WORK_DIR

echo "CP2K calculation completed. Check cp2k.log for results." >> $WORK_DIR/run.log