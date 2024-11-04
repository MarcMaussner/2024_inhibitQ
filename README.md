# InhibitQ Project

This repository contains quantum chemistry calculations and quantum computing simulations for corrosion inhibitor molecules on aluminum surfaces. The repository is organized into two main sections: phase2_submission (containing files for official submission) and calculational_details (containing supporting research and development work that we carried out during the project).

## Repository Structure

### Phase 2 Submission
The `phase2_submission` folder contains the final calculations and results for two inhibitor molecules:
- Al_triazole_qiskit_adaptVQE: Calculations using adaptive VQE approach
- Al_triazole_qiskit_onlyVQE: Calculations using standard VQE approach
- Al_Triazole-3-thiol_qiskit_adaptVQE: Calculations for thiol variant

Each molecule calculation follows a systematic workflow:
1. Classical calculations (folder: 0_classical_calculations)
2. Supercell calculations (folder: 2_supercell)
3. Aluminum surface calculations (folder: 3_Al)
4. Inhibitor molecule calculations (folder: 4_inhibitor)

### Calculational Details
The `calculational_details` folder contains:
- Proof-of-concept implementations
- Research documentation
- Code examples and tutorials
- Resource estimation calculations

## Technical Details

### Software Stack
- CP2K (version 2024.3) for quantum chemistry calculations
- Qiskit for quantum computing simulations
- CUDA Quantum for GPU-accelerated quantum simulations (still in development)

### Running Calculations

#### CP2K Calculations
You can run CP2K calculations using either local installation or Docker:
