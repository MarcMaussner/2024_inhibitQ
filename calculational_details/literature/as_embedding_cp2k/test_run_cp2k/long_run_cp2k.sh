#!/bin/bash

# Ensure script is executed with bash
if [ -z "$BASH_VERSION" ]; then
    echo "This script must be run with bash."
    exit 1
fi

# Set environment variables to avoid shared memory issues
export OMPI_MCA_btl=^vader
export OMPI_MCA_btl=tcp,self
export OMPI_MCA_btl_base_verbose=100
export TMPDIR=/tmp

# Check if MPICH is installed, install if necessary
if ! command -v mpirun &> /dev/null; then
    echo "MPICH not found, installing MPICH..."
    brew install mpich
    brew unlink open-mpi
    brew link mpich
fi

# Verify installation
if ! command -v mpirun &> /dev/null; then
    echo "Failed to install MPICH."
    exit 1
fi

# Ensure CP2K is installed via Homebrew
if ! command -v cp2k.psmp &> /dev/null; then
    echo "CP2K not found, installing CP2K..."
    brew install cp2k
fi

# Verify CP2K installation
if ! command -v cp2k.psmp &> /dev/null; then
    echo "Failed to install CP2K."
    exit 1
fi

# Run CP2K with the provided input file
if [ $# -eq 0 ]; then
    echo "No input file provided. Usage: ./run_cp2k.sh <input_file>"
    exit 1
fi

INPUT_FILE=$1
BASENAME=$(basename "$INPUT_FILE" .inp)
OUTPUT_FILE="${BASENAME}.out"
RESTART_FILES=("${BASENAME}-1.restart" "${BASENAME}-1.restart.bak-1" "${BASENAME}-BFGS.Hessian" "${BASENAME}-pos-1.xyz")

if [ ! -f "$INPUT_FILE" ]; then
    echo "Input file $INPUT_FILE does not exist."
    exit 1
fi

echo "Running CP2K with input file: $INPUT_FILE"
mpirun -np 4 cp2k.psmp -i $INPUT_FILE -o $OUTPUT_FILE

if [ $? -eq 0 ]; then
    echo "CP2K run completed successfully. Output saved to $OUTPUT_FILE"
    # Create a directory with the same name as the input file (without extension) and move the restart files into it
    mkdir -p "$BASENAME"
    mv $OUTPUT_FILE "${RESTART_FILES[@]}" "$BASENAME/"
else
    echo "CP2K run failed. Check the output for details."
    echo "Displaying last 20 lines of the output file for debugging:"
    tail -n 20 $OUTPUT_FILE
fi
