# #run_experiment_cpu.sh
#!/bin/bash
set -e
export OMP_NUM_THREADS=2
export CP2K_DATA_DIR=/opt/cp2k/data
rm *.log *.pdos *.cube || true
echo "Starting CP2K..."
#mpiexec -np 8 cp2k -i MgO.inp > cp2k_output.log 2>&1 &
#mpiexec -np 18 cp2k -i supp.inp > cp2k_output.log 2>&1 &
mpiexec -np 8 cp2k -i Al111.inp > cp2k_output.log 2>&1 &
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