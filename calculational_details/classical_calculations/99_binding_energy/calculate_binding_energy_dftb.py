from ase.io import read
from ase.calculators.dftb import Dftb
import argparse

def setup_dftb_calculator():
    """Setup DFTB+ calculator with D3 correction."""
    calc = Dftb(
        label='binding_energy',
        Hamiltonian_='DFTB',
        Hamiltonian_MaxAngularMomentum_='',
        Hamiltonian_MaxAngularMomentum_Al='p',
        Hamiltonian_MaxAngularMomentum_C='p',
        Hamiltonian_MaxAngularMomentum_H='s',
        Hamiltonian_MaxAngularMomentum_N='p',
        Hamiltonian_Dispersion='DftD3',
        Hamiltonian_DftD3_Damping='BJ',
        Hamiltonian_DftD3_s6='1.0',
        Hamiltonian_DftD3_s8='0.5883',
        Hamiltonian_DftD3_a1='0.5719',
        Hamiltonian_DftD3_a2='3.6017'
    )
    return calc

def calculate_energy(atoms, calculator):
    """Calculate single point energy."""
    atoms.calc = calculator
    return atoms.get_potential_energy()

def main():
    parser = argparse.ArgumentParser(description='Calculate binding energy using DFTB+ with D3 correction.')
    args = parser.parse_args()

    # Read structures
    print("Reading structures...")
    supercell = read('supercell.xyz', format='extxyz')
    substrate = read('substrate.xyz', format='extxyz')
    inhibitor = read('inhibitor.xyz', format='extxyz')

    # Setup calculator
    print("Setting up DFTB+ calculator...")
    calc = setup_dftb_calculator()

    # Calculate energies
    print("Calculating energies...")
    E_supercell = calculate_energy(supercell, calc)
    E_substrate = calculate_energy(substrate, calc)
    E_inhibitor = calculate_energy(inhibitor, calc)
    
    # Calculate binding energy
    E_binding = E_supercell - (E_substrate + E_inhibitor)
    
    # Print results
    print("\nEnergy Results:")
    print(f"Supercell energy:  {E_supercell:.6f} eV")
    print(f"Substrate energy:  {E_substrate:.6f} eV")
    print(f"Inhibitor energy:  {E_inhibitor:.6f} eV")
    print(f"\nBinding Energy:")
    print(f"E_binding = {E_binding:.6f} eV")
    print(f"E_binding = {E_binding/27.211386245988:.6f} Hartree")
    print(f"E_binding = {E_binding * 23.061} kcal/mol")

    # Save results
    with open("binding_energy_results.txt", "w") as f:
        f.write("Energy Results (eV):\n")
        f.write(f"Supercell energy:  {E_supercell:.6f}\n")
        f.write(f"Substrate energy:  {E_substrate:.6f}\n")
        f.write(f"Inhibitor energy:  {E_inhibitor:.6f}\n")
        f.write(f"\nBinding Energy:\n")
        f.write(f"E_binding = {E_binding:.6f} eV\n")
        f.write(f"E_binding = {E_binding/27.211386245988:.6f} Hartree\n")
        f.write(f"E_binding = {E_binding * 23.061} kcal/mol\n")

if __name__ == "__main__":
    main() 