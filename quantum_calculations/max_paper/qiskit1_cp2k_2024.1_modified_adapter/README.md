To run the CPU version:
```bash
docker build -t cp2k-qiskit-cpu .
```
To run the container:
```bash
docker run -d --name cp2k-qiskit-cpu-container -v ${PWD}:/workspace cp2k-qiskit-cpu
```

nohup ./run_experiment_cpu.sh >> run_outcome.log 2>&1 &
pgrep -fl run_experiment_cpu.sh
pgrep -fl cp2k

## Content:
### determine active orbitals:
- to get the active orbitals needed for the active space
### hybrid calculation working:
- working Al111 hybrid calculation but stuck after connection in active space (AS)
### hybrid calculation advanced:
- for step-by-step doing HF and other things as in MgO example in the respective subfolder hybrid-calculation-MgO
