# pylint: disable=invalid-name #test from container
"""CP2K + Qiskit Nature embedding.

Usage:
    python client-vqe-ucc.py
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
from qiskit_algorithms.optimizers import L_BFGS_B
from qiskit.primitives import Estimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer.primitives import Estimator as AerEstimator

from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.algorithms.excited_states_solvers import QEOM, EvaluationRule
from qiskit_nature.second_q.circuit.library import UCC, HartreeFock
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.problems import ElectronicStructureProblem

from qiskit_nature_cp2k.cp2k_integration import CP2KIntegration

# Custom implementation of StatefulVQE and StatefulAdaptVQE
class StatefulVQE:
    def __init__(self, estimator, ansatz, optimizer):
        self.estimator = estimator
        self.ansatz = ansatz
        self.optimizer = optimizer
        self.initial_point = None

class StatefulAdaptVQE:
    def __init__(self, solver, eigenvalue_threshold, gradient_threshold, max_iterations):
        self.solver = solver
        self.eigenvalue_threshold = eigenvalue_threshold
        self.gradient_threshold = gradient_threshold
        self.max_iterations = max_iterations

np.set_printoptions(linewidth=500, precision=6, suppress=True)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

if __name__ == "__main__":
    HOST = "embedding_socket"
    PORT = 12345
    UNIX = False

    parser = argparse.ArgumentParser()
    parser.add_argument("--nalpha", type=int, default=None)
    parser.add_argument("--nbeta", type=int, default=None)
    parser.add_argument("--norbs", type=int, default=None)
    parser.add_argument("--two-qubit-reduce", action="store_true")
    parser.add_argument("--adapt", action="store_true")
    parser.add_argument("--aer", action="store_true")
    args = parser.parse_args()

    if args.nalpha is None or args.nbeta is None or args.norbs is None:
        raise ValueError("Missing argument!")

    num_alpha, num_beta = args.nalpha, args.nbeta
    num_orbs = args.norbs

    if args.two_qubit_reduce:
        mapper = ParityMapper(num_particles=(num_alpha, num_beta))
    else:
        mapper = ParityMapper()

    # Create a dummy problem for initialization
    dummy_problem = ElectronicStructureProblem(num_spatial_orbitals=num_orbs, num_particles=(num_alpha, num_beta))

    initial_state = HartreeFock(
        dummy_problem.num_spatial_orbitals,
        dummy_problem.num_particles,
        mapper,
    )
    ansatz = UCC(
        dummy_problem.num_spatial_orbitals,
        dummy_problem.num_particles,
        "sd",
        mapper,
        initial_state=initial_state,
    )

    def _no_fail(*args, **kwargs):
        return True

    ansatz._check_ucc_configuration = _no_fail

    if args.adapt:
        operator_pool = []
        for op in ansatz.operators:
            for pauli, coeff in zip(op.paulis, op.coeffs):
                if sum(pauli.x & pauli.z) % 2 == 0:
                    continue
                operator_pool.append(SparsePauliOp([pauli], coeffs=[coeff]))

        ansatz = EvolvedOperatorAnsatz(
            operators=operator_pool,
            initial_state=initial_state,
        )

    if args.aer:
        estimator = AerEstimator(approximation=True)
    else:
        estimator = Estimator()

    optimizer = L_BFGS_B()
    solver = StatefulVQE(estimator, ansatz, optimizer)
    solver.initial_point = [0.0] * ansatz.num_parameters

    if args.adapt:
        solver = StatefulAdaptVQE(
            solver,
            eigenvalue_threshold=1e-4,
            gradient_threshold=1e-4,
            max_iterations=1,
        )

    algo = GroundStateEigensolver(mapper, solver)

    integ = CP2KIntegration(algo)
    integ.connect_to_socket(HOST, PORT, UNIX)
    integ.run()
    problem = integ.construct_problem()

    def my_generator(num_spatial_orbitals, num_particles):
        singles = FermionicOp.gen_excitation_list(1, num_spatial_orbitals, num_particles, preserve_spin=False)
        doubles = []
        return singles + doubles

    if isinstance(integ.algo.solver, StatefulAdaptVQE):
        algo = GroundStateEigensolver(integ.algo.qubit_mapper, integ.algo.solver.solver)

    qeom = QEOM(
        algo,
        estimator,
        my_generator,
        aux_eval_rules=EvaluationRule.ALL,
    )

    logger.info(
        "Removing the ElectronicDensity property before starting QEOM"
    )
    problem.properties.electronic_density = None

    excited_state_result = qeom.solve(problem)

    logger.info("QEOM result: \n\n%s\n", excited_state_result)

    logger.info(
        "Excitation Energies:\n\n%s\n",
        excited_state_result.raw_result.excitation_energies,
    )
    logger.info("Transition Amplitudes")
    for (
        key,
        values,
    ) in excited_state_result.raw_result.transition_amplitudes.items():
        logger.info(key)
        for name, val in values.items():
            logger.info(f"\t{name}: {val[0]}")