from ase.build import fcc111
from ase.io import write
from ase.optimize import BFGS
from ase.constraints import FixAtoms

# Assuming MGLCalculator does not depend on DGL directly, otherwise you need to replace it with an alternative
from matgl.ext.ase import MGLCalculator

def create_al_slab(size, layers=4, vacuum=10.0, a=4.05):
    """
    Create an aluminum (111) slab with specified size, layers, and vacuum.
    
    Parameters:
    size (tuple): (x, y) size of the slab, e.g. (2,2)
    layers (int): Number of atomic layers in the slab (default is 4)
    vacuum (float): Vacuum size in Angstroms to add above and below the slab (default is 10.0)
    a (float): Lattice constant for aluminum in Angstroms (default is 4.05)
    
    Returns:
    ase.Atoms: The created Al(111) slab
    """
    slab = fcc111('Al', size=(size[0], size[1], layers), a=a, vacuum=vacuum)
    
    # Save the initial structure
    filename = f"al_slab_{size[0]}x{size[1]}x{layers}_v{vacuum}.xyz"
    write(filename, slab)
    print(f"Created {filename} with {len(slab)} atoms")
    
    return slab

def relax_slab_matgl(slab, fmax=0.05, fix_bottom_layers=2):
    """
    Relax the given slab structure using MatGL potential.
    
    Parameters:
    slab (ase.Atoms): The slab to relax
    fmax (float): The maximum force allowed on any atom for the relaxation to be considered
                  converged (default is 0.05 eV/Angstrom)
    fix_bottom_layers (int): Number of bottom layers to fix during relaxation (default is 2)
    
    Returns:
    ase.Atoms: The relaxed slab structure
    """
    # Set up the MatGL calculator without DGL
    calc = MGLCalculator(model="M3GNet-MP-2021.2.8-PES", device="cpu")
    slab.calc = calc
    
    # Fix the bottom layers
    if fix_bottom_layers > 0:
        mask = [atom.tag > fix_bottom_layers for atom in slab]
        slab.set_constraint(FixAtoms(mask=mask))
    
    # Relax the structure
    dyn = BFGS(slab, trajectory='al_slab_relaxation_matgl.traj')
    dyn.run(fmax=fmax)
    
    # Save the relaxed structure
    relaxed_filename = 'al_slab_relaxed_matgl.xyz'
    write(relaxed_filename, slab)
    print(f"Created relaxed {relaxed_filename} with {len(slab)} atoms")
    
    return slab

def create_and_relax_al_slabs(sizes, layers=4, vacuum=10.0, a=4.05, fmax=0.05):
    """
    Create and relax aluminum (111) slabs of different sizes using MatGL potential.
    
    Parameters:
    sizes (list of tuples): List of (x, y) sizes for the slab, e.g. [(2,2), (3,3), (4,4)]
    layers (int): Number of atomic layers in the slab (default is 4)
    vacuum (float): Vacuum size in Angstroms to add above and below the slab (default is 10.0)
    a (float): Lattice constant for aluminum in Angstroms (default is 4.05)
    fmax (float): The maximum force allowed on any atom for the relaxation to be considered
                  converged (default is 0.05 eV/Angstrom)
    """
    for size in sizes:
        # Create the slab
        slab = create_al_slab(size, layers, vacuum, a)
        
        # Relax the slab with MatGL potential
        relaxed_slab = relax_slab_matgl(slab, fmax)
        
        # Clean up to free memory
        del slab, relaxed_slab

if __name__ == '__main__':
    # Create and relax 2x2x4, 3x3x4, 4x4x4, and 5x5x4 slabs (4 layers each) with 10 Å vacuum
    sizes_to_create = [(2,2), (3,3), (4,4), (5,5)]
    create_and_relax_al_slabs(sizes_to_create, layers=4, vacuum=10.0)
