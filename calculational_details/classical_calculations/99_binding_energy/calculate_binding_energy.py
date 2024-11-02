import ase
from ase.io import read
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
import argparse

def calculate_energy(atoms, calculator):
    """Calculate single point energy for a structure."""
    atoms.calc = calculator
    return atoms.get_potential_energy()

def main():
    parser = argparse.ArgumentParser(description='Calculate binding energy using ORB-D3.')
    parser.add_argument('--device', default='cpu', choices=['cpu', 'cuda'], 
                       help='Device to run calculations on (default: cpu)')
    args = parser.parse_args()

    # Initialize ORB model
    print(f"Initializing ORB-D3 model on {args.device}...")
    orbff = pretrained.orb_d3_v2(device=args.device)
    calc = ORBCalculator(orbff, device=args.device)

    # Read structures
    print("Reading structures...")
    supercell = read('supercell.xyz', format='extxyz')
    substrate = read('substrate.xyz', format='extxyz')
    inhibitor = read('inhibitor.xyz', format='extxyz')

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