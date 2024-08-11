#!/bin/bash
export OMP_NUM_THREADS=1
export CP2K_DATA_DIR=/usr/local/cp2k/data
export CUDA_VISIBLE_DEVICES=0
mpirun -np 4 cp2k.psmp -i MgO.inp > cp2k_output.log 2>&1 &
CP2K_PID=$!
sleep 30
python3 client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5
wait $CP2K_PID
echo "Experiment completed"