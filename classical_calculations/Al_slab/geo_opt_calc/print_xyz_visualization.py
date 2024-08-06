import sys
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase.visualize.plot import plot_atoms

def visualize_xyz(filename):
    """
    Read an XYZ file (standard or extended) and create visualizations.
    
    Parameters:
    filename (str): Path to the XYZ file
    
    Returns:
    None, but saves a PNG visualization of the structure
    """
    # Read the XYZ file
    atoms = read(filename)
    
    # Generate output filename
    png_filename = f"{filename.rsplit('.', 1)[0]}_visualization.png"
    
    # Create visualizations of the structure (top and side views)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    
    # Top view
    plot_atoms(atoms, ax1, radii=0.5, rotation=('0x,0y,0z'))
    ax1.set_title('Top View')
    ax1.axis('off')
    
    # Side view
    plot_atoms(atoms, ax2, radii=0.5, rotation=('90x,0y,0z'))
    ax2.set_title('Side View')
    ax2.axis('off')
    
    plt.tight_layout()
    plt.savefig(png_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved visualization as {png_filename}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python print_xyz_visualization.py <path_to_xyz_file>")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    visualize_xyz(xyz_file)