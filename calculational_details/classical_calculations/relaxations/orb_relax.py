from ase.io import read, write
from ase.optimize import BFGS
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator

def optimize_geometry(input_xyz, output_xyz, fmax=0.005, steps=300):
    # Read the input XYZ file
    atoms = read(input_xyz)

    # Set up the ORB model and calculator
    device = "cpu"  # Change to "cuda" if you have a GPU and want to use it
    orbff = pretrained.orb_d3_v1(device=device)
    calc = ORBCalculator(orbff, device=device)

    # Set the calculator for the atoms
    atoms.set_calculator(calc)

    # Perform geometry optimization
    optimizer = BFGS(atoms)
    optimizer.run(fmax=fmax, steps=steps)

    # Save the optimized structure
    write(output_xyz, atoms)

    print(f"Optimization completed. Final energy: {atoms.get_potential_energy()} eV")
    print(f"Optimized structure saved to {output_xyz}")

# Example usage
input_file = "al_slab_with_triazole_4x4x6_v10.0.xyz"
output_file = "al_slab_with_triazole_4x4x6_v10.0_relaxed.xyz"
optimize_geometry(input_file, output_file)