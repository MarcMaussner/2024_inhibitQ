# Phase 2 Submission - InhibitQ Project

This folder contains the quantum chemistry and quantum computing calculations for the Phase 2 submission of the InhibitQ project, focusing on corrosion inhibitor molecules on aluminum surfaces.

## Directory Structure

The calculations are organized for two inhibitor molecules:
1. Al_triazole_qiskit_adaptVQE (Using adaptive VQE approach)
2. Al_triazole_qiskit_onlyVQE (Using standard VQE approach)
3. Al_Triazole-3-thiol_qiskit_adaptVQE (Thiol variant with adaptive VQE)

Each molecule's calculations follow a systematic workflow with the following subdirectories:
- `0_classical_calculations/`: Initial classical DFT calculations
- `2_supercell/`: Periodic boundary condition calculations
- `3_Al/`: Aluminum surface calculations
- `4_inhibitor/`: Inhibitor molecule calculations

## Technical Details

### Quantum Chemistry Calculations (CP2K)
- Software Version: CP2K 2024.3
- DFT Settings:
  - Spin restricted Kohn-Sham (RKS)
  - PBE functional
  - Periodic boundary conditions
  - Electronic temperature: 1000K
  - Cutoff settings:
    - Density: 1.0E-10
    - Gradient: 1.0E-10
    - Tau: 1.0E-10

### Quantum Computing Simulations
- Framework: Qiskit
- VQE Implementations:
  - Standard VQE
  - Adaptive VQE (with dynamic operator pool)
- Output files:
  - `cp2k.log`: CP2K calculation results
  - `python_output.log`: VQE calculation results

## Running the Calculations

Each calculation directory contains a `run.sh` script that executes both CP2K and VQE calculations, example usage can be:
```bash
cd Al_triazole_qiskit_adaptVQE/2_supercell
nohup ./run.sh &
```

## References

1. CP2K: [J. Chem. Phys. 152, 194103 (2020)](https://doi.org/10.1063/5.0007045)
2. Quickstep: [Comput. Phys. Commun. 167, 103-128 (2005)](https://doi.org/10.1016/j.cpc.2004.12.014)
3. PBE Functional: Perdew, Burke, Ernzerhof, Phys. Rev. Letter, vol. 77, n 18, pp. 3865-3868 (1996)

## Notes
- All calculations use restricted Kohn-Sham DFT
- The VQE calculations are performed after the DFT calculations complete
- Molecular geometries and electronic structure data are preserved in respective output files
