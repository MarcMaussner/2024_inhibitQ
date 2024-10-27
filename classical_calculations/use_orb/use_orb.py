import ase
from ase.io import read, write
import numpy as np
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
from ase.optimize import BFGS
from ase import Atoms
import re
import sys
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Optimize structure using ORB.')
    parser.add_argument('input_file', help='Path to the input XYZ file')
    parser.add_argument('--device', default='cpu', choices=['cpu', 'cuda'], 
                       help='Device to run calculations on (default: cpu)')
    parser.add_argument('--fmax', type=float, default=0.01,
                       help='Maximum force criterion for optimization (default: 0.01)')
    parser.add_argument('--output', default='optimized_structure.xyz',
                       help='Output file name (default: optimized_structure.xyz)')
    
    args = parser.parse_args()

    # Parse the extended XYZ format
    def parse_lattice(lattice_str):
        # Extract numbers from the lattice string
        numbers = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", lattice_str)]
        # Convert to 3x3 matrix
        return np.array(numbers).reshape(3, 3)

    # Read the structure from the file
    try:
        with open(args.input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Could not find file '{args.input_file}'")
        sys.exit(1)

    # Parse header
    try:
        header = lines[1]
        lattice_str = re.search(r'Lattice="([^"]+)"', header).group(1)
        pbc_str = re.search(r'pbc="([^"]+)"', header).group(1)
    except (IndexError, AttributeError):
        print("Error: File format incorrect. Expected extended XYZ format with Lattice and pbc information.")
        sys.exit(1)

    # Convert lattice string to cell matrix
    cell = parse_lattice(lattice_str)

    # Convert pbc string to boolean array
    pbc = [x.strip().lower() == 't' for x in pbc_str.split()]

    # Parse atomic positions and create lists for atoms
    symbols = []
    positions = []
    tags = []

    for line in lines[2:]:
        if line.strip():
            try:
                parts = line.split()
                symbols.append(parts[0])
                positions.append([float(x) for x in parts[1:4]])
                tags.append(int(parts[4]))
            except (IndexError, ValueError) as e:
                print(f"Error parsing line: {line.strip()}")
                print(f"Error details: {e}")
                sys.exit(1)

    # Create ASE Atoms object
    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=pbc,
        tags=tags
    )

    print(f"Successfully loaded structure with {len(atoms)} atoms")

    # Initialize the ORB model and calculator
    print(f"Initializing ORB model on {args.device}...")
    orbff = pretrained.orb_d3_v2(device=args.device)
    calc = ORBCalculator(orbff, device=args.device)
    atoms.set_calculator(calc)

    # Print initial energy
    initial_energy = atoms.get_potential_energy()
    print(f"Initial energy: {initial_energy:.3f} eV")

    # Run geometry optimization
    print("Starting geometry optimization...")
    opt = BFGS(atoms, trajectory='optimization.traj')
    opt.run(fmax=args.fmax)

    # Print final energy
    final_energy = atoms.get_potential_energy()
    print(f"Final energy: {final_energy:.3f} eV")
    print(f"Energy change: {final_energy - initial_energy:.3f} eV")

    # Save optimized structure
    write(args.output, atoms, format='extxyz')
    print(f"Optimized structure saved to {args.output}")

if __name__ == "__main__":
    main()