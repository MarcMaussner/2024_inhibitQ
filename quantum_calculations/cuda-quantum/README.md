To download CUDA Quantum, you can use Docker with this command:

```sh
docker pull nvcr.io/nvidia/nightly/cuda-quantum:latest
```

There is an example for running Water Molecule with Active Space (CPU vs. GPU) available at this [link](https://nvidia.github.io/cuda-quantum/latest/examples/python/tutorials/vqe_water_active_space.html).

Then you can run the image to generate a container like this

```sh
docker run --gpus all -it --name cuda-quantum20240801 -v $(pwd):/workspace nvcr.io/nvidia/nightly/cuda-quantum:latest
```

To monitor the nvidia GPU execution, you can run
```sh
watch -d -n 1 nvidia-smi
```

For cuQuantum, you can do the following:
```sh
docker pull nvcr.io/nvidia/cuquantum-appliance:24.03-cuda12.2.2-devel-ubuntu22.04-x86_64
docker run --gpus all -it --name cuQuantum24.03 -v $(pwd):/workspace nvcr.io/nvidia/cuquantum-appliance:24.03-cuda12.2.2-devel-ubuntu22.04-x86_64
```

