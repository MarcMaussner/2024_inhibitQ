# InhibitQ Project

**Authors:**
- Hazem Abdelhafiz [hazem.AbdElhafiz@gmail.com](mailto:hazem.AbdElhafiz@gmail.com)
- Karim Elgammal [egkarim@gmail.com](mailto:egkarim@gmail.com)
- Marc Maußner [marc.maussner@infoteam.de](mailto:marc.maussner@infoteam.de)

*Authors listed in alphabetical order.*

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

## How to Cite

### Dataset Citation
If you use this dataset in your research, please cite:
```bibtex
@misc{inhibitq2024dataset,
title={InhibitQ: Quantum Computing Simulations of Corrosion Inhibitors},
author={Abdelhafiz, Hazem and Elgammal, Karim and Maußner, Marc},
year={2024},
doi={10.24435/materialscloud:e2-dr},
url={https://doi.org/10.24435/materialscloud:e2-dr}
}
```


### Associated Publication
This dataset supports the following publication:
- Karim Elgammal and Marc Maußner (2024). A Quantum Computing Approach to Simulating Corrosion Inhibition. *arXiv preprint* arXiv:2412.00951.
  - Preprint: https://arxiv.org/abs/2412.00951
  - Field: Condensed Matter > Materials Science

```bibtex
@misc{elgammal2024,
title={A Quantum Computing Approach to Simulating Corrosion Inhibition},
author={Karim Elgammal and Marc Maußner},
year={2024},
eprint={2412.00951},
archivePrefix={arXiv},
primaryClass={cond-mat.mtrl-sci},
url={https://arxiv.org/abs/2412.00951},
}

```