from ase.build import fcc111
from ase.io import write
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

def create_al_slabs(sizes, layers=4, vacuum=10.0, a=4.05):
    """
    Create aluminum (111) slabs of different sizes with multiple layers and vacuum.
    
    Parameters:
    sizes (list of tuples): List of (x, y) sizes for the slab, e.g. [(2,2), (3,3), (4,4)]
    layers (int): Number of atomic layers in the slab (default is 4)
    vacuum (float): Vacuum size in Angstroms to add above and below the slab (default is 10.0)
    a (float): Lattice constant for aluminum in Angstroms (default is 4.05)
    
    Returns:
    None, but saves .xyz files and .png visualizations for each slab
    """
    for size in sizes:
        # Create the slab with specified number of layers and vacuum
        slab = fcc111('Al', size=(size[0], size[1], layers), a=a, vacuum=vacuum)
        
        # Generate filenames
        xyz_filename = f"al_slab_{size[0]}x{size[1]}x{layers}_v{vacuum}.xyz"
        png_filename = f"al_slab_{size[0]}x{size[1]}x{layers}_v{vacuum}.png"
        
        # Save the slab as an XYZ file
        write(xyz_filename, slab)
        
        # Create visualizations of the slab (top and side views)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # Top view
        plot_atoms(slab, ax1, radii=0.5, rotation=('0x,0y,0z'))
        ax1.set_title('Top View')
        ax1.axis('off')
        
        # Side view
        plot_atoms(slab, ax2, radii=0.5, rotation=('90x,0y,0z'))
        ax2.set_title('Side View')
        ax2.axis('off')
        
        plt.tight_layout()
        plt.savefig(png_filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created {xyz_filename} with {len(slab)} atoms")
        print(f"Saved visualization as {png_filename}")

# Example usage
sizes_to_create = [(4,4)]  #[(2,2), (3,3)]
create_al_slabs(sizes_to_create, layers=6, vacuum=10.0)