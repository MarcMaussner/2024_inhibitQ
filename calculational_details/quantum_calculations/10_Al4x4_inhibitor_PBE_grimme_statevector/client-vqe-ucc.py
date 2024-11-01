# pylint: disable=invalid-name
"""CP2K + Qiskit Nature embedding.

Usage:
    python dft-emb-client.py
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
from scipy.sparse import issparse
from scipy.sparse.linalg import eigsh, eigs
from qiskit.primitives import Estimator
from qiskit_nature.logging import logging as nature_logging
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.algorithms.excited_states_solvers import (
    QEOM, EvaluationRule)
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit_nature.second_q.circuit.library.ansatzes.utils import \
    generate_fermionic_excitations
from qiskit_nature.second_q.mappers import ParityMapper

from qiskit_nature_cp2k.cp2k_integration import CP2KIntegration
from qiskit_algorithms import NumPyEigensolver, NumPyMinimumEigensolver

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

# Custom NumPyEigensolver to handle sparse matrices
class CustomNumPyEigensolver(NumPyEigensolver):
    def _solve_sparse(self, op_matrix, k):
        n = op_matrix.shape[0]
        if k is None:
            k = n
        if self._filter_criterion is not None:
            # If we have a filter criterion, we need to solve for all eigenvalues
            k = n

        # Check if the matrix is Hermitian
        is_hermitian = np.allclose(op_matrix.todense(), op_matrix.todense().conj().T)

        if is_hermitian:
            eigvals, eigvecs = eigsh(op_matrix, k=min(k, n-1), which='SA')
        else:
            eigvals, eigvecs = eigs(op_matrix, k=min(k, n-1), which='SR')

        # Sort eigenvalues and eigenvectors
        idx = eigvals.argsort()
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

        return eigvals, eigvecs

class CustomNumPyMinimumEigensolver(NumPyMinimumEigensolver):
    def __init__(self):
        super().__init__()
        self._eigensolver = CustomNumPyEigensolver()

if __name__ == "__main__":
    HOST = "embedding_socket"
    PORT = 12345
    UNIX = True

    parser = argparse.ArgumentParser()
    parser.add_argument("--nalpha", type=int, default=None)
    parser.add_argument("--nbeta", type=int, default=None)
    parser.add_argument("--norbs", type=int, default=None)
    parser.add_argument("--two-qubit-reduce", action="store_true")
    args = parser.parse_args()

    if args.nalpha is None or args.nbeta is None or args.norbs is None:
        raise ValueError("Missing argument!")

    num_alpha, num_beta = args.nalpha, args.nbeta
    num_orbs = args.norbs

    if args.two_qubit_reduce:
        mapper = ParityMapper(num_particles=(num_alpha, num_beta))
    else:
        mapper = ParityMapper()

    initial_state = HartreeFock(
        num_orbs,
        (num_alpha, num_beta),
        mapper,
    )

    numpy_solver = CustomNumPyMinimumEigensolver()
    algo = GroundStateEigensolver(mapper, numpy_solver)

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

    estimator = Estimator()

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