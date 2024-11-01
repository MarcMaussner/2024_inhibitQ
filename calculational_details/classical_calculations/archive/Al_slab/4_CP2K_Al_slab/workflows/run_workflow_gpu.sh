#!/bin/bash

export CP2K_TAG=v2023.2

# Run CP2K calculations
docker run -it --rm --gpus all -v $PWD:/workspace -w /workspace nvcr.io/hpc/cp2k:${CP2K_TAG} python3 cp2k_calc.py

# Run analysis
python3 analysis.py