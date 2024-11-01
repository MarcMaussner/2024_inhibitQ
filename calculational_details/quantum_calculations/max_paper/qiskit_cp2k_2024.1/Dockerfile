FROM cp2k/cp2k:2024.1_openmpi_generic_psmp

# install python and OpenMPI
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git \
    openmpi-bin \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir --upgrade pip

# Install qiskit-nature-cp2k after changing from 'qiskit.algorithms' to 'qiskit_algorithms'
RUN git clone -b fix-qiskit1.x-compatibility https://github.com/KarimElgammal/qiskit-nature-cp2k.git && \
    cd qiskit-nature-cp2k && \
    pip3 install .

RUN pip3 install qiskit-aer mpi4py

# Set up working directory
WORKDIR /workspace

# Copy necessary files
COPY MgO.inp /workspace/
COPY MgO.xyz /workspace/
COPY client-vqe-ucc.py /workspace/
COPY run_calculation.sh /workspace/

# Make sure the run script is executable
RUN chmod +x /workspace/run_calculation.sh

# Set environment variables
ENV CP2K_DATA_DIR=/opt/cp2k/data
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Keep the container running
CMD ["tail", "-f", "/dev/null"]