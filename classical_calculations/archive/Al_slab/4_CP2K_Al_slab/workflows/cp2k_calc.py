import os
from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.io import read, write
from ase.calculators.cp2k import CP2K
from ase.optimize import BFGS
from ase.constraints import FixAtoms

def create_slab(size=(2,2,1), vacuum=20, a=4.05):
    """Create an Al(111) slab."""
    slab = fcc111('Al', size=size, a=a, vacuum=vacuum)
    return slab

def add_molecule(slab, molecule_file, height=2.0):
    """Add a molecule to the slab."""
    molecule = read(molecule_file)
    add_adsorbate(slab, molecule, height=height, position=(0, 0))
    return slab

def setup_cp2k_calc(label='cp2k_calc', cutoff=600, rel_cutoff=60, basis_set='DZVP-MOLOPT-SR-GTH', 
                    pseudo_potential='GTH-PBE'):
    """Set up a CP2K calculator with vdW-DF2 functional."""
    inp = f"""
    &FORCE_EVAL
      &DFT
        &XC
          &XC_FUNCTIONAL
            &VDW_POTENTIAL
              POTENTIAL_TYPE NON_LOCAL
              &NON_LOCAL
                TYPE VDW_DF2
              &END NON_LOCAL
            &END VDW_POTENTIAL
          &END XC_FUNCTIONAL
        &END XC
        &PRINT
          &PDOS
            NLUMO 10
            COMPONENTS
          &END PDOS
        &END PRINT
      &END DFT
      &PRINT
        &FORCES ON
      &END PRINT
    &END FORCE_EVAL
    &GLOBAL
      PRINT_LEVEL MEDIUM
    &END GLOBAL
    """
    
    calc = CP2K(label=label,
                command='mpirun cp2k_shell',
                basis_set=basis_set,
                pseudo_potential=pseudo_potential,
                max_scf=200,
                cutoff=cutoff,
                rel_cutoff=rel_cutoff,
                inp=inp)
    return calc

def relax_structure(atoms, fmax=0.01, steps=200, trajectory='relaxation.traj'):
    """Relax the given structure using CP2K with vdW-DF2 functional."""
    # Fix bottom two layers
    constraint = FixAtoms(mask=[atom.tag <= 2 for atom in atoms])
    atoms.set_constraint(constraint)

    calc = setup_cp2k_calc()
    atoms.calc = calc

    dyn = BFGS(atoms, trajectory=trajectory)
    dyn.run(fmax=fmax, steps=steps)

    return atoms

def single_point_calculation(atoms):
    """Perform a single point calculation."""
    calc = setup_cp2k_calc()
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    return energy

def main():
    # Create slab
    slab = create_slab(size=(2,2,1), vacuum=20)
    write('slab.xyz', slab)

    # Add adsorbate (if needed)
    # Uncomment the following lines if you want to add an adsorbate
    # adsorbate_file = 'adsorbate.xyz'
    # slab_with_adsorbate = add_molecule(slab, adsorbate_file)
    # write('slab_with_adsorbate.xyz', slab_with_adsorbate)

    # Relax structure
    relaxed_slab = relax_structure(slab)
    write('relaxed_slab.xyz', relaxed_slab)

    # Perform single point calculation
    energy = single_point_calculation(relaxed_slab)
    print(f"Total energy: {energy} eV")

if __name__ == "__main__":
    main()