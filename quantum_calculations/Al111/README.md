use command:
bash
docker run -v $PWD:/mnt --shm-size=1g -it --rm cp2k/cp2k mpiexec -genv OMP_NUM_THREADS=2 -np 18 cp2k Al111.inp
