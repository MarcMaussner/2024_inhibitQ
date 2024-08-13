# Running 2024.2 version
```sh
docker build --no-cache -t cp2k-qiskit-cpu . && 
docker run -d --name cp2k-qiskit-cpu-container -v ${PWD}:/workspace cp2k-qiskit-cpu
```