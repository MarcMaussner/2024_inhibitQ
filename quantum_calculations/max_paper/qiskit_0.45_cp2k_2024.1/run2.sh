#!/bin/bash
set -e

NUM_MPI_PROCESSES=48
export OMP_NUM_THREADS=2
export OMP_STACKSIZE=1G

rm -f *.log
rm -f embedding_socket

echo "Starting CP2K with MPI..."

mpirun -np ${NUM_MPI_PROCESSES} cp2k -i MgO.inp > cp2k_output.log 2>&1 &
CP2K_PID=$!

echo "CP2K PID: $CP2K_PID"
echo "Waiting for CP2K to initialize..."

socket_created=false
for i in {1..300}; do
    if [ -S embedding_socket ]; then
        echo "Socket file created after $i seconds"
        socket_created=true
        break
    fi
    if [ $((i % 60)) -eq 0 ]; then
        echo "CP2K process status after $i seconds:"
        ps -p $CP2K_PID -o pid,cmd,state,start,etime,time
        echo "Last 20 lines of CP2K output:"
        tail -n 20 cp2k_output.log
    fi
    sleep 1
done

if ! $socket_created; then
    echo "Socket file not created after 300 seconds. CP2K may still be running."
    echo "Last 50 lines of CP2K output:"
    tail -n 50 cp2k_output.log
    echo "Continuing with the Python client..."
fi

echo "Starting the VQE client..."
python3 -u client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5 2>&1 | tee python_output.log

echo "Waiting for CP2K to finish..."
wait $CP2K_PID

echo "Experiment completed"
echo "CP2K output:"
cat cp2k_output.log

echo "Python client output:"
cat python_output.log