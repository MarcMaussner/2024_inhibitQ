# Use an image with Python 3.12
FROM python:3.12-slim

# Set environment variables to prevent Python from writing pyc files to disc and buffering stdout and stderr
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    git \
    cmake \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge (a minimal conda installer for ARM architecture)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O Miniforge3-Linux-aarch64.sh \
    && bash Miniforge3-Linux-aarch64.sh -b -p /opt/conda \
    && rm Miniforge3-Linux-aarch64.sh

# Add conda to PATH
ENV PATH=/opt/conda/bin:$PATH

# Install ASE, DFTD4, and CHGNet using conda and pip
RUN conda install -c conda-forge ase \
    && conda install -c conda-forge dftd4 \
    && pip install chgnet

# Set the working directory
WORKDIR /usr/src/app

# Copy the entrypoint script and Python script
COPY entrypoint.sh .
COPY chgnet_vdw_script_corrected.py .

# Make the entrypoint script executable
RUN chmod +x entrypoint.sh

# Define the entry point
ENTRYPOINT ["./entrypoint.sh"]
