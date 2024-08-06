# 2024_inhibitQ
- repository to store papers, thoughts and proof-of-concept examples for the 2024 Airbus Coating challenge

## Repository Organization

The repository is organized as follows:

- `papers/`: Contains research papers and related documents.
- `thoughts/`: A collection of notes, ideas, and brainstorming sessions.
- `proof-of-concept/`: Contains proof-of-concept examples and code.
  - `example1/`: Description of example 1.
  - `example2/`: Description of example 2.
- `scripts/`: Utility scripts for various tasks.
- `data/`: Datasets used in the project.
- `results/`: Results from experiments and simulations.

## Running CP2K using NVIDIA NGC Docker

To run CP2K using the NVIDIA NGC Docker container, you can use the following command:

```sh
export CP2K_TAG=v2023.2
docker run -it --rm --gpus all -v ${PWD}:/workspace -w /workspace nvcr.io/hpc/cp2k:${CP2K_TAG} cp2k.psmp -i H2O-dft-ls.NREP2.inp

You can update your [`README.md`](command:_github.copilot.openRelativePath?%5B%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2FUsers%2Fkarim%2Fgithub%2F2024_bmw_airbus%2FREADME.md%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%5D "/Users/karim/github/2024_bmw_airbus/README.md") file to include information about how the repository is organized. Here's an example of how you can structure it:

```markdown
# 2024_inhibitQ
- repository to store papers, thoughts and proof-of-concept examples for the 2024 Airbus Coating challenge

## Repository Organization

The repository is organized as follows:

- `papers/`: Contains research papers and related documents.
- `thoughts/`: A collection of notes, ideas, and brainstorming sessions.
- `proof-of-concept/`: Contains proof-of-concept examples and code.
  - `example1/`: Description of example 1.
  - `example2/`: Description of example 2.
- `scripts/`: Utility scripts for various tasks.
- `data/`: Datasets used in the project.
- `results/`: Results from experiments and simulations.
- `literature/`: Contains a walkthrough on how to run CP2K locally.

## Running Local CP2K Calculations

To run local CP2K calculations using the `run_cp2k.sh` script, follow these steps:

1. Ensure you have the necessary dependencies installed:
   - Bash
   - MPICH
   - CP2K

2. If MPICH is not installed, the script will attempt to install it using Homebrew.

3. If CP2K is not installed, the script will attempt to install it using Homebrew.

4. Run the script with the input file as an argument:
   ```sh
   ./run_cp2k.sh <input_file>
   ```
Replace <input_file> with the path to your CP2K input file.

The script will execute CP2K with the provided input file and save the output to a file with the same base name but with a .out extension.

If the run is successful, the output and restart files will be moved to a directory named after the input file (without the extension).

Example:

This will run CP2K with the H2O-dft-ls.NREP2.inp input file and save the output to H2O-dft-ls.NREP2.out.

## Running CP2K using NVIDIA NGC Docker

To run CP2K using the NVIDIA NGC Docker container, you can use the following command:

```sh
export CP2K_TAG=v2023.2
docker run -it --rm --gpus all -v ${PWD}:/workspace -w /workspace nvcr.io/hpc/cp2k:${CP2K_TAG} cp2k.psmp -i H2O-dft-ls.NREP2.inp
```

Example input files can be downloaded from [this link](https://github.com/cp2k/cp2k/blob/c415b5ddc864fe89e1e8e74ebbc33ad8b898175d/benchmarks/QS_DM_LS/README.md).

## Downloading CUDA Quantum

To download CUDA Quantum, you can use Docker with this command:

```sh
docker pull nvcr.io/nvidia/nightly/cuda-quantum:latest
```

There is an example for running Water Molecule with Active Space (CPU vs. GPU) available at [this link](https://nvidia.github.io/cuda-quantum/latest/examples/python/tutorials/vqe_water_active_space.html).
```

note: package to visualize the results of the simulations
https://www.ovito.org/