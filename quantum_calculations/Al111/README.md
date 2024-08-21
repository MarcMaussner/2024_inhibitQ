To run the Al111 use command:

Latest working commands:
```bash
python -u client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5 2>&1 | tee python_output.log
docker run -v $PWD:/mnt --shm-size=1g -u $(id -u):$(id -g) -it --rm cp2k/cp2k mpiexec -genv OMP_NUM_THREADS=2 -np 18 cp2k Al111_active_space.inp
```