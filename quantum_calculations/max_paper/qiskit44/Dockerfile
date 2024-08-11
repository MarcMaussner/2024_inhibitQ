FROM nvcr.io/hpc/cp2k:v2023.2

# Install Python, pip, and git
RUN apt-get update && apt-get install -y python3 python3-pip git

# Install specific versions of Qiskit and other required packages
RUN pip3 install qiskit==0.44.0 qiskit-aer==0.12.0 qiskit-nature==0.6.0

# Install specific version of qiskit-algorithms from GitHub
RUN git clone -b 0.3.0 https://github.com/qiskit-community/qiskit-algorithms.git && \
    cd qiskit-algorithms && \
    pip3 install .

# Clone and install your forked qiskit-nature-cp2k
RUN git clone -b fix-package-metadata https://github.com/KarimElgammal/qiskit-nature-cp2k.git && \
    cd qiskit-nature-cp2k && \
    pip3 install .

# Set working directory
WORKDIR /workspace

# Copy your scripts and input files
COPY MgO.inp ./
COPY MgO.xyz ./
COPY client-vqe-ucc.py ./
COPY run_experiment.sh ./

# Create keep_alive.sh in the container
RUN echo '#!/bin/bash\ntail -f /dev/null' > /keep_alive.sh && \
    chmod +x /keep_alive.sh

# Make sure the scripts are executable
RUN chmod +x run_experiment.sh

# Set the entry point
ENTRYPOINT ["/keep_alive.sh"]