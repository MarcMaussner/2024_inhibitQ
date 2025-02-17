# #run_experiment_cpu.sh
#!/bin/bash
set -e
export OMP_NUM_THREADS=2
export CP2K_DATA_DIR=/opt/cp2k/data
rm *.log *.pdos *.cube embedding_socket || true
echo "Starting CP2K..."

mpiexec -np 12 cp2k -i MgO.inp > cp2k_output.log 2>&1 &
CP2K_PID=$!

echo "Waiting for CP2K to initialize..."
sleep 30

echo "Starting Python client..."
python3 -u client-vqe-ucc.py --nalpha 3 --nbeta 3 --norbs 5 2>&1 | tee python_output.log

echo "Waiting for CP2K to finish..."
wait $CP2K_PID

echo "Experiment completed"
echo "CP2K output:"
cat cp2k_output.log
echo "Python output:"
cat python_output.log