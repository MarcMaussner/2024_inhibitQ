#!/usr/bin/env python3
# pylint: disable=invalid-name
"""CP2K + Qiskit Nature embedding.

Usage:
    python client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5
"""

import time
import os

def wait_for_socket(socket_path, timeout=300, retry_interval=5):
    start_time = time.time()
    while time.time() - start_time < timeout:
        if os.path.exists(socket_path):
            return True
        time.sleep(retry_interval)
    return False

import sys

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

logger.info("Python script started")
logger.info(f"Python version: {sys.version}")
logger.info(f"Current working directory: {os.getcwd()}")
logger.info(f"Contents of current directory: {os.listdir()}")

from __future__ import annotations

import argparse
import logging
import os
import time
import socket

import numpy as np
from qiskit.algorithms.optimizers import L_BFGS_B
from qiskit.circuit.library import EvolvedOperatorAnsatz
from qiskit.primitives import Estimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer.primitives import Estimator as AerEstimator
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.algorithms.excited_states_solvers import QEOM, EvaluationRule
from qiskit_nature.second_q.circuit.library import UCC, HartreeFock
from qiskit_nature.second_q.operators.fermionic_op import FermionicOp
from qiskit_nature.second_q.mappers import ParityMapper

from qiskit_nature_cp2k.cp2k_integration import CP2KIntegration
from qiskit_nature_cp2k.stateful_adapt_vqe import StatefulAdaptVQE
from qiskit_nature_cp2k.stateful_vqe import StatefulVQE

np.set_printoptions(linewidth=500, precision=6, suppress=True)

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def wait_for_socket(socket_path, timeout=300):
    start_time = time.time()
    while not os.path.exists(socket_path):
        if time.time() - start_time > timeout:
            logger.error(f"Socket {socket_path} not created after {timeout} seconds")
            return False
        time.sleep(1)
    logger.info(f"Socket {socket_path} found after {time.time() - start_time:.2f} seconds")
    return True

def my_generator(num_spatial_orbitals, num_particles):
    singles = FermionicOp.generate_fermionic_excitations(
        1, num_spatial_orbitals, num_particles, generalized=True
    )
    return singles


logger.info("Entering main execution block")
if __name__ == "__main__":
    HOST = "embedding_socket"
    PORT = 12345
    UNIX = True

    parser = argparse.ArgumentParser()
    parser.add_argument("--nalpha", type=int, required=True)
    parser.add_argument("--nbeta", type=int, required=True)
    parser.add_argument("--norbs", type=int, required=True)
    parser.add_argument("--two-qubit-reduce", action="store_true")
    parser.add_argument("--adapt", action="store_true")
    parser.add_argument("--aer", action="store_true")
    args = parser.parse_args()

    num_alpha, num_beta = args.nalpha, args.nbeta
    num_orbs = args.norbs

    logger.info(f"Starting calculation with nalpha={num_alpha}, nbeta={num_beta}, norbs={num_orbs}")
    logger.debug(f"Arguments: {args}")

    mapper = ParityMapper(num_particles=(num_alpha, num_beta)) if args.two_qubit_reduce else ParityMapper()

    initial_state = HartreeFock(num_orbs, (num_alpha, num_beta), mapper)
    ansatz = UCC(
        num_orbs,
        (num_alpha, num_beta),
        "sd",
        mapper,
        initial_state=initial_state,
    )

    if args.adapt:
        logger.info("Using ADAPT-VQE")
        operator_pool = [
            SparsePauliOp([pauli], coeffs=[coeff])
            for op in ansatz.operators
            for pauli, coeff in zip(op.paulis, op.coeffs)
            if sum(pauli.x & pauli.z) % 2 != 0
        ]
        ansatz = EvolvedOperatorAnsatz(operators=operator_pool, initial_state=initial_state)

    estimator = AerEstimator(approximation=True) if args.aer else Estimator()

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

    logger.debug(f"Current working directory: {os.getcwd()}")
    logger.debug(f"Contents of current directory: {os.listdir()}")

    socket_path = os.path.join(os.getcwd(), HOST)
    logger.debug(f"Looking for socket at: {socket_path}")

    if not wait_for_socket(socket_path):
        logger.error("Socket not found. Exiting.")
        exit(1)

    logger.info("Attempting to connect to socket...")
    integ = CP2KIntegration(algo)
    try:
        integ.connect_to_socket(HOST, PORT, UNIX)
        logger.info("Successfully connected to socket")
    except Exception as e:
        logger.exception(f"Failed to connect to socket: {e}")
        # Try to get more information about the socket
        try:
            sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            sock.connect(socket_path)
            logger.info("Successfully created a new connection to the socket")
            sock.close()
        except Exception as sock_e:
            logger.exception(f"Failed to create a new connection to the socket: {sock_e}")
        exit(1)

    logger.info("Running CP2K integration...")
    try:
        integ.run()
        logger.info("CP2K integration completed successfully")
    except Exception as e:
        logger.exception(f"Error during CP2K integration: {e}")
        exit(1)

    logger.info("Constructing problem...")
    try:
        problem = integ.construct_problem()
        logger.info("Problem constructed successfully")
    except Exception as e:
        logger.exception(f"Error constructing problem: {e}")
        exit(1)

    logger.info("Solving with QEOM...")
    if isinstance(integ.algo.solver, StatefulAdaptVQE):
        algo = GroundStateEigensolver(integ.algo.qubit_mapper, integ.algo.solver.solver)

    qeom = QEOM(
        algo,
        estimator,
        my_generator,
        aux_eval_rules=EvaluationRule.ALL,
    )

    logger.info("Removing the ElectronicDensity property before starting QEOM")
    problem.properties.electronic_density = None

    try:
        excited_state_result = qeom.solve(problem)
        logger.info("QEOM calculation completed successfully")
    except Exception as e:
        logger.exception(f"Error during QEOM calculation: {e}")
        exit(1)

    logger.info("QEOM result: \n\n%s\n", excited_state_result)

    logger.info(
        "Excitation Energies:\n\n%s\n",
        excited_state_result.raw_result.excitation_energies,
    )
    logger.info("Transition Amplitudes")
    for key, values in excited_state_result.raw_result.transition_amplitudes.items():
        logger.info(key)
        for name, val in values.items():
            logger.info(f"\t{name}: {val[0]}")

    logger.info("Calculation completed successfully")