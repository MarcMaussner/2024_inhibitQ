Right now, the calculation is converging at the active space calculation
PBE functional used + Grimme

Steps the client-vqe-ucc.py does:
1) Create parser for arguments
2) Create ParityMapper
3) Set initial_state to HartreeFock
4) Create Ansatz UCC
5) If AdaptVQE is used:
   - Create EvolvedOperatorAnsatz (with reduced PauliOps)
6) Create estimator:
   either AerEstimator, Estimator or StatevectorEstimator
7) Create solver (StatefulVQE)
8) If AdaptVQE -> create solver (StatefulAdaptVQE)
9) IF Numpy -> create solver (NumPyMinimumEigensolver)
10) Start CP2K integration
11) Create algo:
    either local or cp2k GroundStateEigensolver
12) create QEOM with fermionic_excitations (see: https://arxiv.org/abs/1910.12890)
13) solve problem with QEOM
14) Print results