# pylint: disable=invalid-name
"""CP2K + Qiskit Nature embedding.

Usage:
    python dft-emb-client.py
"""

from __future__ import annotations

import argparse
import json
import logging

import numpy as np
from scipy.sparse import isspmatrix
from scipy.sparse.linalg import eigsh
from qiskit.primitives import Estimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer.primitives import Estimator as AerEstimator
from qiskit_nature.logging import logging as nature_logging
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.algorithms.excited_states_solvers import (
    QEOM, EvaluationRule)
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit_nature.second_q.circuit.library.ansatzes.utils import \
    generate_fermionic_excitations
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_algorithms.eigensolvers.numpy_eigensolver import NumPyEigensolver

from qiskit_nature_cp2k.cp2k_integration import CP2KIntegration

np.set_printoptions(linewidth=500, precision=6, suppress=True)

logger = logging.getLogger(__name__)

level = logging.DEBUG

nature_logging.set_levels_for_names(
    {
        __name__: level,
        "qiskit": level,
        "qiskit_nature": level,
        "qiskit_nature_cp2k": level,
    }
)

class CustomNumPyEigensolver(NumPyEigensolver):
    @staticmethod
    def _solve_sparse(op_matrix, k=None):
        # Check if matrix is Hermitian
        is_hermitian = (op_matrix != op_matrix.conj().T).nnz == 0

        if is_hermitian:
            # Use eigsh for Hermitian matrices
            if k is None or k >= op_matrix.shape[0] - 1:
                k = op_matrix.shape[0] - 1
            eigval, eigvec = eigsh(op_matrix, k=k, which="SA")
        else:
            # For non-Hermitian, fall back to dense eigenvalue solver
            dense_matrix = op_matrix.todense()
            eigval, eigvec = np.linalg.eig(dense_matrix)
            
        if k is not None:
            eigval = eigval[:k]
            eigvec = eigvec[:, :k]
        
        return eigval, eigvec

if __name__ == "__main__":
    HOST = "embedding_socket"
    PORT = 12345
    UNIX = True

    parser = argparse.ArgumentParser()
    parser.add_argument("--nalpha", type=int, default=None)
    parser.add_argument("--nbeta", type=int, default=None)
    parser.add_argument("--norbs", type=int, default=None)
    args = parser.parse_args()

    if args.nalpha is None or args.nbeta is None or args.norbs is None:
        raise ValueError("Missing argument!")

    num_alpha, num_beta = args.nalpha, args.nbeta
    num_orbs = args.norbs

    mapper = ParityMapper()

    initial_state = HartreeFock(
        num_orbs,
        (num_alpha, num_beta),
        mapper,
    )

    solver = CustomNumPyEigensolver()
    algo = GroundStateEigensolver(mapper, solver)

    integ = CP2KIntegration(algo)
    integ.connect_to_socket(HOST, PORT, UNIX)
    integ.run()
    problem = integ.construct_problem()

    def my_generator(num_spatial_orbitals, num_particles):
        singles = generate_fermionic_excitations(
            1, num_spatial_orbitals, num_particles, preserve_spin=False
        )
        doubles = []
        return singles + doubles

    estimator = AerEstimator(run_options={"shots": None})  # None shots means using statevector

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