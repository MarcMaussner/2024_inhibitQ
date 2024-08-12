#!/bin/bash
set -e
export OMP_NUM_THREADS=1
export CP2K_DATA_DIR=/usr/local/cp2k/data
export CUDA_VISIBLE_DEVICES=0

echo "Starting CP2K..."
mpirun -np 4 cp2k.psmp -i MgO.inp > cp2k_output.log 2>&1 &
CP2K_PID=$!

echo "Waiting for CP2K to initialize..."
sleep 30

echo "Starting Python client..."
python3 -u client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5

echo "Waiting for CP2K to finish..."
wait $CP2K_PID

echo "Experiment completed"
echo "CP2K output:"
cat cp2k_output.log