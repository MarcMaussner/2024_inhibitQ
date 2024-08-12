FROM nvcr.io/hpc/cp2k:v2023.2

# Install Python, pip, and git
RUN apt-get update && apt-get install -y python3 python3-pip git

# Install Qiskit and other required packages
RUN pip3 install qiskit qiskit-aer qiskit-nature

# Install latest qiskit-algorithms from GitHub
RUN pip3 install git+https://github.com/qiskit-community/qiskit-algorithms.git
# RUN git clone https://github.com/qiskit-community/qiskit-algorithms.git && \
#     cd qiskit-algorithms && \
#     pip3 install .

# Clone and install your forked qiskit-nature-cp2k
RUN git clone -b fix-qiskit1.x-compatibility https://github.com/KarimElgammal/qiskit-nature-cp2k.git && \
    cd qiskit-nature-cp2k && \
    pip3 install .

# Set working directory
WORKDIR /workspace

# Create a script to keep the container running
RUN echo '#!/bin/bash\ntail -f /dev/null' > /keep_alive.sh && chmod +x /keep_alive.sh

# Set the entry point to the keep-alive script
ENTRYPOINT ["/keep_alive.sh"]

# Copy your scripts and input files (commented out)
# COPY MgO.inp ./
# COPY MgO.xyz ./
# COPY client-vqe-ucc.py ./
# COPY run_experiment.sh ./

