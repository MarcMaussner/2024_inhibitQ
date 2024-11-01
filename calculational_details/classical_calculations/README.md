## Running CP2K using NVIDIA NGC Docker

To run CP2K using the NVIDIA NGC Docker container, you can use the following command:

```sh
export CP2K_TAG=v2023.2
docker run -it --rm --gpus all -v ${PWD}:/workspace -w /workspace nvcr.io/hpc/cp2k:${CP2K_TAG} cp2k.psmp -i H2O-dft-ls.NREP2.inp

Example input files can be downloaded from this [link](https://github.com/cp2k/cp2k/blob/c415b5ddc864fe89e1e8e74ebbc33ad8b898175d/benchmarks/QS_DM_LS/README.md).
```

To run cp2k on CPU, you can use the following commands:
```sh
docker run -it --rm cp2k/cp2k:latest
```

or with OpenMPI with:
```sh
docker run -it --rm --gpus all --shm-size=1g -v $PWD:/mnt -u $(id -u $USER):$(id -g $USER) --name cp2k_2024.1 cp2k/cp2k:2024.1_openmpi_generic_psmp mpirun -bind-to none -np 4 -x OMP_NUM_THREADS=2 cp2k -i H2O-dft-ls.NREP2.inp
```

You can check more details in this [link](https://github.com/cp2k/cp2k-containers?tab=readme-ov-file)