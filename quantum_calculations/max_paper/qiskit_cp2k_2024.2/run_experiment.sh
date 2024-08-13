#!/bin/bash
set -e

# Set the number of MPI processes and OpenMP threads
NUM_MPI_PROCESSES= 4 #6
export OMP_NUM_THREADS= 2 #16

# Ensure CP2K can find its data files
export CP2K_DATA_DIR=/opt/cp2k/data

# Clean up any existing log files
rm *.log || true

echo "Starting CP2K..."
mpirun -np ${NUM_MPI_PROCESSES} cp2k -i MgO.inp > cp2k_output.log 2>&1 &
CP2K_PID=$!

echo "Waiting for CP2K to initialize..."
sleep 30

echo "Starting Python client..."
python3 -u client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5 2>&1 | tee python_output.log

echo "Waiting for CP2K to finish..."
wait $CP2K_PID

echo "Experiment completed"
echo "CP2K output:"
cat cp2k_output.log
echo "Python output:"
cat python_output.log