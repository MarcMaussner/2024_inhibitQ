#!/bin/bash
# run_docker.sh

# Build the Docker image
docker build -t chgnet_dftd4 .

# Run the Docker container, mounting the current directory to /usr/src/app in the container
docker run --rm -v $(pwd):/usr/src/app chgnet_dftd4
