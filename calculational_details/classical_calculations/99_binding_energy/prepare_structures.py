import ase
from ase.io import read, write
from ase import Atoms
import numpy as np
import re
import sys

def read_extended_xyz(filename):
    """Read extended XYZ file with proper handling of cell and PBC."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    header = lines[1]
    lattice_str = re.search(r'Lattice="([^"]+)"', header).group(1)
    pbc_str = re.search(r'pbc="([^"]+)"', header).group(1)
    
    numbers = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", lattice_str)]
    cell = np.array(numbers).reshape(3, 3)
    pbc = [x.strip().lower() == 't' for x in pbc_str.split()]
    
    symbols = []
    positions = []
    tags = []
    
    for line in lines[2:]:
        if line.strip():
            parts = line.split()
            symbols.append(parts[0])
            positions.append([float(x) for x in parts[1:4]])
            tags.append(int(parts[4]))
    
    return Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=pbc,
        tags=tags
    )

def main():
    # Read the supercell structure
    supercell = read_extended_xyz("al_slab_with_triazole_4x4x6_v10.0.xyz")
    
    # Create substrate (keep the cell parameters)
    substrate_indices = [i for i, tag in enumerate(supercell.get_tags()) if tag > 0]
    substrate = Atoms(
        symbols=[supercell.symbols[i] for i in substrate_indices],
        positions=[supercell.positions[i] for i in substrate_indices],
        cell=supercell.get_cell(),
        pbc=supercell.get_pbc(),
        tags=[supercell.get_tags()[i] for i in substrate_indices]
    )
    
    # Create inhibitor (keep the cell parameters)
    inhibitor_indices = [i for i, tag in enumerate(supercell.get_tags()) if tag == 0]
    inhibitor = Atoms(
        symbols=[supercell.symbols[i] for i in inhibitor_indices],
        positions=[supercell.positions[i] for i in inhibitor_indices],
        cell=supercell.get_cell(),
        pbc=supercell.get_pbc(),
        tags=[supercell.get_tags()[i] for i in inhibitor_indices]
    )
    
    # Save structures
    write('substrate.xyz', substrate, format='extxyz')
    write('inhibitor.xyz', inhibitor, format='extxyz')
    write('supercell.xyz', supercell, format='extxyz')
    
    # Print structure information
    print("\nStructure Information:")
    print(f"Supercell: {len(supercell)} atoms")
    print(f"Cell parameters: {supercell.cell.diagonal()}")
    print(f"PBC: {supercell.get_pbc()}")
    print(f"\nSubstrate: {len(substrate)} atoms")
    print(f"Elements: {set(substrate.get_chemical_symbols())}")
    print(f"\nInhibitor: {len(inhibitor)} atoms")
    print(f"Elements: {set(inhibitor.get_chemical_symbols())}")

if __name__ == "__main__":
    main() 