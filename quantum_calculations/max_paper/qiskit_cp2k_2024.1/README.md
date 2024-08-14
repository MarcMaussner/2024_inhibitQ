# Running 2024.1 version
```sh
docker build -t cp2k2024.1-qiskit .
docker run -d --name cp2k2024.1-qiskit-container -v $(pwd):/workspace cp2k2024.1-qiskit
```