# 2024_inhibitQ
- repository to store papers, thoughts and proof-of-concept examples for the 2024 Airbus Coating challenge

## the python environment
ubuntu@ip-172-31-36-110:~/qiskit-nature-cp2k$ pip list
Package                           Version     Editable project location
--------------------------------- ----------- -------------------------------
amazon-braket-default-simulator   1.26.0
amazon-braket-schemas             1.22.1
amazon-braket-sdk                 1.88.1
annotated-types                   0.7.0
antlr4-python3-runtime            4.9.2
backoff                           2.2.1
backports.entry-points-selectable 1.3.0
boltons                           24.0.0
boto3                             1.35.49
botocore                          1.35.49
certifi                           2024.8.30
charset-normalizer                3.4.0
cloudpickle                       3.1.0
decorator                         5.1.1
dill                              0.3.9
h5py                              3.12.1
idna                              3.10
importlib_metadata                8.5.0
jmespath                          1.0.1
mpmath                            1.3.0
mypy-extensions                   1.0.0
nest-asyncio                      1.6.0
networkx                          3.4.2
numpy                             2.1.2
openpulse                         1.0.1
openqasm3                         1.0.0
opt_einsum                        3.4.0
oqpy                              0.3.7
pbr                               6.1.0
pip                               24.3.1
psutil                            6.1.0
py                                1.11.0
pydantic                          2.9.2
pydantic_core                     2.23.4
python-dateutil                   2.9.0.post0
python-dotenv                     1.0.1
qiskit                            1.2.4
qiskit-aer                        0.15.1
qiskit-algorithms                 0.3.1
qiskit-braket-provider            0.4.1
qiskit-ionq                       0.5.8
qiskit-nature                     0.7.2
qiskit_nature_cp2k                0.0.1       /home/ubuntu/qiskit-nature-cp2k
requests                          2.32.3
retry                             0.9.2
rustworkx                         0.15.1
s3transfer                        0.10.3
scipy                             1.14.1
setuptools                        75.2.0
six                               1.16.0
stevedore                         5.3.0
symengine                         0.13.0
sympy                             1.13.3
typing_extensions                 4.12.2
urllib3                           2.2.3
zipp                              3.20.2
ubuntu@ip-172-31-36-110:~/qiskit-nature-cp2k$ 

# OLD text

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

testing