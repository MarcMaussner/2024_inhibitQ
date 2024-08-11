# pylint: disable=invalid-name
"""CP2K + Qiskit Nature embedding.

Usage:
    python client-vqe-ucc.py
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
from qiskit.algorithms.optimizers import L_BFGS_B
from qiskit.circuit.library import EvolvedOperatorAnsatz
from qiskit.primitives import Estimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer.primitives import Estimator as AerEstimator

from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.algorithms.excited_states_solvers import QEOM, EvaluationRule
from qiskit_nature.second_q.circuit.library import UCC, HartreeFock
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp

from qiskit_nature_cp2k.cp2k_integration import CP2KIntegration
from qiskit_nature_cp2k.stateful_adapt_vqe import StatefulAdaptVQE
from qiskit_nature_cp2k.stateful_vqe import StatefulVQE

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

    mapper = JordanWignerMapper()

    initial_state = HartreeFock(
        num_spin_orbitals=2*num_orbs,
        num_particles=(num_alpha, num_beta),
        qubit_mapper=mapper,
    )
    ansatz = UCC(
        qubit_mapper=mapper,
        num_particles=(num_alpha, num_beta),
        num_spin_orbitals=2*num_orbs,
        excitations='sd',
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