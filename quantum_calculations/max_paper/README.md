This record contains some inputs, outputs and scripts to reproduce the results of the manuscript "A general framework for active space embedding methods: applications in quantum computing".

The files here are regarding some calculations related to VQE
The python script contain the workflow/input for the Qiskit process 
The python scripts provide the utilities for and set up the fermion-to-qubit mapping, generating the quantum circuits and executing them.

- MgO.inp: the input file
- MgO.out: the output file
- MgO.xyz: the (initial) geometry
- MgO-RESTART.wfn: the wave function file to restar the calculation with a good guess
- MgO-kx-y.pdos: the projected density of states
- MgO-pos-1.xyz: all the geometries at intermediate optimization steps
- qiskit.log: the output from Qiskit
- MgO-1_0.fcidump: one- and two-electron integrals over the active orbitals defining the embedding space


To run the container and be able to attach it to vs code
```sh
docker build --no-cache -t cp2k-qiskit . #you can remove no-cache for faster execution
docker build -t cp2k-qiskit . # for faster build
docker run -d --name cp2k-qiskit-container --gpus all -v ${PWD}:/workspace cp2k-qiskit ./run_experiment.sh
docker run -d --name cp2k-qiskit-container --gpus all -v ${PWD}:/workspace cp2k-qiskit #if you want to attach it to vs code
```

