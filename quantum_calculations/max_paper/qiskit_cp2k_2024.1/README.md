# Running 2024.1 version
```sh
docker build -t cp2k-qiskit-2024.1 .

docker run -it --name cp2k-qiskit-2024.1-container --shm-size=32g \
  --cpuset-cpus="0-95" \
  -v ${PWD}:/workspace \
  cp2k-qiskit-2024.1 /bin/bash
  ```