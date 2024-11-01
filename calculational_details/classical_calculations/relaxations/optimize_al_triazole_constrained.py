from ase.io import read, write
from ase.optimize import BFGS
from ase.constraints import FixAtoms, FixedPlane
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
import numpy as np

def optimize_geometry(input_xyz, output_xyz, fmax=0.01, steps=200, max_height=25.0):
    # Read the input XYZ file
    atoms = read(input_xyz)

    # Ensure PBC is set correctly
    atoms.set_pbc([True, True, False])

    # Fix the bottom four layers of the aluminum slab
    z_positions = atoms.positions[:, 2]
    unique_z = np.sort(np.unique(z_positions))
    if len(unique_z) > 4:
        z_cutoff = unique_z[3]  # Fix atoms below the fourth unique z-position
        fixed_indices = [atom.index for atom in atoms if atom.position[2] < z_cutoff]
        constraint_fixed = FixAtoms(indices=fixed_indices)

    # Add FixedPlane constraint to keep the molecule close to the surface
    slab_top = unique_z[3]  # Top of the slab (4th layer)
    molecule_indices = [atom.index for atom in atoms if atom.position[2] > slab_top]
    constraints = [constraint_fixed]
    for idx in molecule_indices:
        constraint_plane = FixedPlane(idx, [0, 0, 1])  # Constrain motion in z-direction
        constraints.append(constraint_plane)

    atoms.set_constraint(constraints)

    # Set up the ORB model and calculator
    device = "cpu"  # Change to "cuda" if you have a GPU and want to use it
    orbff = pretrained.orb_d3_v1(device=device)
    calc = ORBCalculator(orbff, device=device)

    # Set the calculator for the atoms
    atoms.set_calculator(calc)

    # Perform geometry optimization
    optimizer = BFGS(atoms, trajectory='optimization.traj')
    
    def check_height(atoms=atoms):
        max_z = np.max(atoms.positions[:, 2])
        return max_z <= max_height

    optimizer.attach(check_height, interval=1)
    optimizer.run(fmax=fmax, steps=steps)

    # Save the optimized structure
    write(output_xyz, atoms)

    print(f"Optimization completed. Final energy: {atoms.get_potential_energy()} eV")
    print(f"Optimized structure saved to {output_xyz}")

# Example usage
input_file = "al_slab_with_triazole_4x4x6_v10.0.xyz"
output_file = "optimized_structure.xyz"
optimize_geometry(input_file, output_file, max_height=26.0)