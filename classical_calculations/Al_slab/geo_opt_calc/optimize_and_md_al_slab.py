import numpy as np
import torch
from ase.io import read, write
from ase.constraints import FixAtoms
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from chgnet.model import CHGNet, StructOptimizer
from chgnet.model.dynamics import MolecularDynamics
from ase.io.trajectory import Trajectory
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_and_prepare_slab(filename, fix_bottom_layers=2):
    """Load the slab from an XYZ file and prepare it for optimization."""
    slab = read(filename)
    
    # Fix bottom layers
    if fix_bottom_layers > 0:
        mask = [atom.index < fix_bottom_layers * len(slab) // 4 for atom in slab]
        slab.set_constraint(FixAtoms(mask=mask))
    
    # Convert to pymatgen Structure
    structure = AseAtomsAdaptor.get_structure(slab)
    
    return structure

def optimize_slab(structure, fmax=0.05):
    """Optimize the slab structure using CHGNet."""
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    logging.info(f"Using device: {device}")

    chgnet = CHGNet.load()
    relaxer = StructOptimizer(model=chgnet, use_device=device)

    try:
        result = relaxer.relax(structure, fmax=fmax)
        optimized_structure = result["final_structure"]
        logging.info("Optimization completed successfully.")
    except Exception as e:
        logging.error(f"Optimization failed: {str(e)}")
        return None

    return optimized_structure

def run_molecular_dynamics(structure, temperature=1000, timestep=2, steps=50):
    """Run molecular dynamics simulation on the optimized structure."""
    chgnet = CHGNet.load()
    
    md = MolecularDynamics(
        atoms=structure,
        model=chgnet,
        ensemble="nvt",
        temperature=temperature,  # in K
        timestep=timestep,  # in femto-seconds
        trajectory="md_out.traj",
        logfile="md_out.log",
        loginterval=100,
    )
    md.run(steps)
    logging.info(f"Completed {steps} steps of MD simulation.")

def save_structure(structure, base_filename):
    """Save structure in CIF, standard XYZ, and extended XYZ formats."""
    # Save as CIF
    structure.to(f"{base_filename}.cif")
    logging.info(f"Saved structure to {base_filename}.cif (CIF)")

    # Convert to ASE Atoms object
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    
    # Save as standard XYZ
    write(f"{base_filename}_standard.xyz", ase_atoms, format='xyz')
    logging.info(f"Saved structure to {base_filename}_standard.xyz (standard XYZ)")
    
    # Save as extended XYZ
    write(f"{base_filename}_extended.xyz", ase_atoms, format='extxyz')
    logging.info(f"Saved structure to {base_filename}_extended.xyz (extended XYZ)")

def main():
    # Load initial slab
    initial_structure = load_and_prepare_slab('Al111.xyz')
    logging.info(f"Loaded initial structure with {len(initial_structure)} atoms.")

    # Optimize slab
    optimized_structure = optimize_slab(initial_structure)

    if optimized_structure is not None:
        # Save optimized structure
        save_structure(optimized_structure, 'Al111_optimized')

        # Run molecular dynamics
        run_molecular_dynamics(optimized_structure)

        # Save final MD structure
        traj = Trajectory("md_out.traj")
        final_md_structure = AseAtomsAdaptor.get_structure(traj[-1])
        save_structure(final_md_structure, 'Al111_after_md')
    else:
        logging.warning("Optimization failed, skipping MD simulation.")

if __name__ == "__main__":
    main()