import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

def analyze_output(filename):
    energy = None
    with open(filename, 'r') as f:
        for line in f:
            if "ENERGY| Total FORCE_EVAL" in line:
                energy = float(line.split()[-1])
    return energy

def analyze_pdos(filename):
    data = np.loadtxt(filename, skiprows=1)
    energy = data[:, 0]
    dos = data[:, 1]
    return energy, dos

def plot_pdos(k1_energy, k1_dos, k2_energy, k2_dos):
    plt.figure(figsize=(10, 6))
    plt.plot(k1_energy, k1_dos, label='k1')
    plt.plot(k2_energy, k2_dos, label='k2')
    plt.xlabel('Energy (eV)')
    plt.ylabel('DOS')
    plt.title('Projected Density of States')
    plt.legend()
    plt.savefig('pdos_plot.png')
    plt.close()

def analyze_xyz(filename):
    atoms = read(filename)
    return atoms

def main():
    # Analyze main output
    if os.path.exists('MgO.out'):
        energy = analyze_output('MgO.out')
        print(f"Final energy for MgO: {energy} a.u.")
    else:
        print("Output file MgO.out not found.")

    # Analyze PDOS files
    if os.path.exists('MgO-k1-1.pdos') and os.path.exists('MgO-k2-1.pdos'):
        k1_energy, k1_dos = analyze_pdos('MgO-k1-1.pdos')
        k2_energy, k2_dos = analyze_pdos('MgO-k2-1.pdos')
        plot_pdos(k1_energy, k1_dos, k2_energy, k2_dos)
        print("PDOS plot saved as pdos_plot.png")
    else:
        print("PDOS files not found.")

    # Analyze XYZ files
    if os.path.exists('MgO-pos-1.xyz'):
        atoms = analyze_xyz('MgO-pos-1.xyz')
        print(f"Number of atoms in final structure: {len(atoms)}")
        print(f"Cell dimensions: {atoms.cell.cellpar()}")
    else:
        print("XYZ file MgO-pos-1.xyz not found.")

    # Check for RESTART file
    if os.path.exists('MgO-RESTART.wfn'):
        print("RESTART file found: MgO-RESTART.wfn")
    else:
        print("RESTART file not found.")

if __name__ == "__main__":
    main()