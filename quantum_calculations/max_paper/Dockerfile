FROM nvcr.io/hpc/cp2k:v2023.2

# Install Python and pip
RUN apt-get update && apt-get install -y python3 python3-pip

# Install Qiskit and other required packages
RUN pip3 install qiskit qiskit-nature

# Clone and install qiskit-nature-cp2k
RUN git clone https://github.com/mrossinek/qiskit-nature-cp2k.git && \
    cd qiskit-nature-cp2k && \
    pip3 install .

# Set working directory
WORKDIR /workspace

# Copy your scripts and input files
COPY MgO.inp /workspace/
COPY MgO.xyz /workspace/
COPY MgO-RESTART.wfn /workspace/
COPY client-vqe-ucc.py /workspace/

# Create a script to run both CP2K and Qiskit
RUN echo '#!/bin/bash \n\
export OMP_NUM_THREADS=1 \n\
export CP2K_DATA_DIR=/usr/local/cp2k/data \n\
export CUDA_VISIBLE_DEVICES=0 \n\
mpirun -np 4 cp2k.psmp -i MgO.inp > cp2k_output.log 2>&1 & \n\
CP2K_PID=$! \n\
sleep 30 \n\
python3 client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5 \n\
wait $CP2K_PID' > /workspace/run_all.sh && chmod +x /workspace/run_all.sh

# Set the entry point
ENTRYPOINT ["/workspace/run_all.sh"]