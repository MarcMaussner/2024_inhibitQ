To run the CPU version:
```bash
docker build -t cp2k-qiskit-cpu .
```
To run the container:
```bash
docker run -d --name cp2k-qiskit-cpu-container -v ${PWD}:/workspace cp2k-qiskit-cpu
```