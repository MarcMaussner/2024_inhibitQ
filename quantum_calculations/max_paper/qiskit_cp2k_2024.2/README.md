# Running 2024.2 version
```sh
docker build -t cp2k-qiskit:latest .
docker run -it --name cp2k-qiskit-container -v ${PWD}:/workspace cp2k-qiskit:latest /bin/bash
docker run -it --name cp2k-qiskit-container --shm-size=2g -v ${PWD}:/workspace cp2k-qiskit:latest /bin/bash

docker build --no-cache -t cp2k-qiskit-cpu . && 
docker run -d --name cp2k-qiskit-cpu-container -v ${PWD}:/workspace cp2k-qiskit-cpu
```

For hpc6a instance:
```sh
docker run -it --name cp2k-qiskit-container --shm-size=32g \
  --cpuset-cpus="0-95" \
  -v ${PWD}:/workspace \
  cp2k-qiskit:latest /bin/bash
  ```