# Contents of this folder is all the preliminary stuff with braket

# Contents:
- VQE_Example:
  - execution of VQE Example for H2 from https://github.com/amazon-braket/amazon-braket-examples/blob/main/examples/hybrid_quantum_algorithms/VQE_Chemistry/VQE_chemistry_braket.ipynb
- first_test:
  - example of runnign local and aws/SV1 bell from https://aws.amazon.com/de/blogs/quantum-computing/setting-up-your-local-development-environment-in-amazon-braket/
- hybrid_hob:
  - example for hybrid job from https://docs.aws.amazon.com/braket/latest/developerguide/braket-hybrid-job-decorator.html
  - hybrid_job_vqe take from https://qiskit-community.github.io/qiskit-braket-provider/tutorials/2_tutorial_hybrid_jobs.html
- track_costs:
  - example for ising model vqe and cost calculation (runnign local simulator) taken from: https://github.com/amazon-braket/amazon-braket-examples/blob/main/examples/hybrid_quantum_algorithms/VQE_Transverse_Ising/VQE_Transverse_Ising_Model.ipynb
- estimate_resources:
  - Notebooks created for resource esimtation of our workflow for requesting further AWS credits
    - old estimates: old files for estimating resources without use of hybrid jobs
    files in folder are using hybridjob feature
- qeom_poc:
  - Notebooks to see if workflow (VQE + qeom for H2) can by applied to braket localSim, SV1, IonQ Aria1, IQM Garnet
