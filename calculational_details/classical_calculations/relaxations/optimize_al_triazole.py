from ase.io import read, write
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
import numpy as np

def optimize_geometry(input_xyz, output_xyz, fmax=0.01, steps=200):
    # Read the input XYZ file
    atoms = read(input_xyz)

    # The cell and PBC information should already be set, but let's ensure it
    if not np.any(atoms.pbc):
        atoms.set_pbc([True, True, False])  # Set PBC to True for x and y, False for z

    # Fix the bottom four layers of the aluminum slab
    z_positions = atoms.positions[:, 2]
    unique_z = np.unique(z_positions)
    if len(unique_z) > 4:
        z_cutoff = unique_z[4]  # Fix atoms below the fifth unique z-position
        constraint = FixAtoms(indices=[atom.index for atom in atoms if atom.position[2] < z_cutoff])
        atoms.set_constraint(constraint)

    # Set up the ORB model and calculator
    device = "cpu"  # Change to "cuda" if you have a GPU and want to use it
    orbff = pretrained.orb_d3_v1(device=device)
    calc = ORBCalculator(orbff, device=device)

    # Set the calculator for the atoms
    atoms.set_calculator(calc)

    # Perform geometry optimization
    optimizer = BFGS(atoms, trajectory='optimization.traj')
    optimizer.run(fmax=fmax, steps=steps)

    # Save the optimized structure
    write(output_xyz, atoms)

    print(f"Optimization completed. Final energy: {atoms.get_potential_energy()} eV")
    print(f"Optimized structure saved to {output_xyz}")

# Example usage
input_file = "al_slab_with_triazole_4x4x6_v10.0.xyz"
output_file = "optimized_structure.xyz"
optimize_geometry(input_file, output_file)