from ase.io import read
from ase.calculators.cp2k import CP2K
import argparse

def setup_cp2k_calculator():
    """Setup CP2K calculator with DFT-D3 correction."""
    calc = CP2K(
        label='binding_energy',
        command='cp2k_shell',
        basis_set='DZVP-MOLOPT-SR-GTH',
        pseudo_potential='GTH-PBE',
        xc='PBE',
        force_eval_method='Quickstep',
        potential_file='GTH_POTENTIALS',
        basis_set_file='BASIS_MOLOPT',
        max_scf=200,
        cutoff=400,
        stress_tensor=True,
        inp="""&DFT
                &VDW_POTENTIAL
                    POTENTIAL_TYPE PAIR_POTENTIAL
                    &PAIR_POTENTIAL
                        TYPE DFTD3
                        PARAMETER_FILE_NAME dftd3.dat
                        REFERENCE_FUNCTIONAL PBE
                    &END PAIR_POTENTIAL
                &END VDW_POTENTIAL
            &END DFT"""
    )
    return calc

def calculate_energy(atoms, calculator):
    """Calculate single point energy."""
    atoms.calc = calculator
    return atoms.get_potential_energy()

def main():
    parser = argparse.ArgumentParser(description='Calculate binding energy using CP2K with DFT-D3.')
    args = parser.parse_args()

    # Read structures
    print("Reading structures...")
    supercell = read('supercell.xyz', format='extxyz')
    substrate = read('substrate.xyz', format='extxyz')
    inhibitor = read('inhibitor.xyz', format='extxyz')

    # Setup calculator
    print("Setting up CP2K calculator...")
    calc = setup_cp2k_calculator()

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