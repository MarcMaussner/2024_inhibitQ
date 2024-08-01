#!/bin/bash

# Ensure the input file is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

input_file=$1
base_name=$(basename "$input_file" .inp)
output_folder="${base_name}_output"

# Create a new directory named after the input file
mkdir -p "$output_folder"

# Copy the input file to the new directory
cp "$input_file" "$output_folder/"

# Change to the new directory
cd "$output_folder"

# Run CP2K with a single process
mpirun -np 1 cp2k.psmp -i "$input_file" -o "${base_name}.out"

# Change back to the original directory
cd ..

echo "CP2K run completed successfully. Output saved to $output_folder"
