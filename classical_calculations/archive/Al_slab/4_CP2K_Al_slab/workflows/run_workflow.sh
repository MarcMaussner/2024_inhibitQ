#!/bin/bash

# Run CP2K calculations
docker run -v $PWD:/workspace -w /workspace cp2k/cp2k python3 cp2k_calc.py

# Run analysis
python3 analysis.py