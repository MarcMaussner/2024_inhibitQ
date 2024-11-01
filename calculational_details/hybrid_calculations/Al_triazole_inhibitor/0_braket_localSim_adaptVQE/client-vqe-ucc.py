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
import socket
import time
import sys

import numpy as np
from qiskit_algorithms.optimizers import L_BFGS_B, SPSA
from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit.circuit.library import EvolvedOperatorAnsatz
from qiskit.primitives import Estimator
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import Aer
from qiskit_aer.primitives import Estimator as AerEstimator
from qiskit_nature.logging import logging as nature_logging
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.algorithms.excited_states_solvers import (
    QEOM, EvaluationRule)
from qiskit_nature.second_q.circuit.library import UCC, HartreeFock
from qiskit_nature.second_q.circuit.library.ansatzes.utils import \
    generate_fermionic_excitations
from qiskit_nature.second_q.mappers import ParityMapper

from qiskit_nature_cp2k.cp2k_integration import CP2KIntegration
from qiskit_nature_cp2k.stateful_adapt_vqe import StatefulAdaptVQE
from qiskit_nature_cp2k.stateful_vqe import StatefulVQE

from braket.aws import AwsDevice
from braket.devices import Devices
from braket.jobs import hybrid_job, save_job_result
from qiskit.primitives import BackendEstimatorV2 as BackendEstimator

from qiskit_braket_provider import BraketProvider
from qiskit_braket_provider import BraketLocalBackend

np.set_printoptions(linewidth=500, precision=6, suppress=True)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('quantum_calculation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

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


if __name__ == "__main__":
    HOST = "embedding_socket"
    PORT = 12345
    UNIX = True

    parser = argparse.ArgumentParser()
    parser.add_argument("--nalpha", type=int, default=None)
    parser.add_argument("--nbeta", type=int, default=None)
    parser.add_argument("--norbs", type=int, default=None)
    parser.add_argument("--two-qubit-reduce", action="store_true")
    parser.add_argument("--adapt", action="store_true")
    parser.add_argument("--aer", action="store_true")
    parser.add_argument("--sv_estimator", action="store_true")
    parser.add_argument("--braket", action="store_true")
    parser.add_argument("--numpy", action="store_true")
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
    ansatz = UCC(
        num_orbs,
        (num_alpha, num_beta),
        "sd",
        mapper,
        # generalized=True,
        # preserve_spin=False,
        initial_state=initial_state,
    )

    def _no_fail(*args, **kwargs):
        return True

    ansatz._check_ucc_configuration = _no_fail

    if args.adapt:
        logger.info("=== ADAPT-VQE mode activated ===")
        operator_pool = []
        for op in ansatz.operators:
            for pauli, coeff in zip(op.paulis, op.coeffs):
                if sum(pauli.x & pauli.z) % 2 == 0:
                    continue
                operator_pool.append(SparsePauliOp([pauli], coeffs=[coeff]))
        
        logger.info(f"Created operator pool with {len(operator_pool)} operators")
        ansatz = EvolvedOperatorAnsatz(
            operators=operator_pool,
            initial_state=initial_state,
        )

    if args.aer:
        # Configure the Aer simulator with the desired number of shots
        backend = Aer.get_backend('aer_simulator')
        estimator = AerEstimator(backend=backend, approximation=True)
    else:
        if args.sv_estimator:
            estimator = StatevectorEstimator()
        else:
            estimator = Estimator()

    if args.braket:
        logger.info("=== Amazon Braket backend activated ===")
        backend = BraketLocalBackend()
        logger.info("Using Braket Local Simulator backend")
        from qiskit.primitives import Estimator as EstimatorV1
        estimator = EstimatorV1()
        estimator.options.default_shots = 1000
        logger.info(f"Configured estimator with {estimator.options.default_shots} shots")

    def callback(nfev, parameters, energy, stepsize):
        logger.info(f"Iteration {nfev}: Energy = {energy:.6f}, Parameters = {parameters}")
        return False

    # Initialize optimizer with more conservative settings
    optimizer = SPSA(
        maxiter=100,
        learning_rate=0.005,
        perturbation=0.05,
        last_avg=1
    )
    logger.info("Configured SPSA optimizer with conservative parameters")

    # Ensure initial point matches the number of parameters
    num_parameters = ansatz.num_parameters
    initial_point = np.zeros(num_parameters)  # Start with zeros instead of random
    logger.info(f"Initialized {num_parameters} parameters to zero")

    solver = StatefulVQE(
        estimator,
        ansatz,
        optimizer,
        callback=callback,
        initial_point=initial_point
    )
    logger.info("Initialized StatefulVQE with matched parameter count")

    if args.adapt:
        class CustomAdaptVQE(StatefulAdaptVQE):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.previous_energies = []
                self.convergence_window = 3
                self.convergence_threshold = 1e-6
                self._theta = np.array([])  # Initialize empty parameter array
                logger.info(f"""ADAPT-VQE Configuration:
                    - Convergence window: {self.convergence_window}
                    - Convergence threshold: {self.convergence_threshold}
                    - Initial parameters: {self._theta}""")

            def solve(self, problem):
                logger.info("Starting ADAPT-VQE solve iteration")
                try:
                    result = super().solve(problem)
                    self.previous_energies.append(result.total_energies[-1])
                    if len(self.previous_energies) > self.convergence_window:
                        self.previous_energies.pop(0)
                    
                    if len(self.previous_energies) == self.convergence_window:
                        energy_range = max(self.previous_energies) - min(self.previous_energies)
                        logger.info(f"Current energy range in window: {energy_range}")
                        if energy_range < self.convergence_threshold:
                            logger.info("ADAPT-VQE Convergence achieved!")
                            self._converged = True
                    
                    return result
                except Exception as e:
                    logger.error(f"Error in ADAPT-VQE solve: {str(e)}")
                    raise

        logger.info("Initializing CustomAdaptVQE solver")
        solver = CustomAdaptVQE(
            solver,
            eigenvalue_threshold=1e-6,
            gradient_threshold=1e-4,
            max_iterations=20,
        )
        logger.info("""ADAPT-VQE solver parameters:
            - Eigenvalue threshold: 1e-6
            - Gradient threshold: 1e-4
            - Maximum iterations: 20""")

    if args.numpy:
        solver = NumPyMinimumEigensolver()

    algo = GroundStateEigensolver(mapper, solver)

    logger.info(
        "Starting CP2KIntegration"
        )
    integ = CP2KIntegration(algo)
    integ.connect_to_socket(HOST, PORT, UNIX)
    integ.run()
    problem = integ.construct_problem()

    def my_generator(num_spatial_orbitals, num_particles):
        singles = generate_fermionic_excitations(
            1, num_spatial_orbitals, num_particles, preserve_spin=False
        )
        doubles = []
        # doubles = generate_fermionic_excitations(
        #     2, num_spatial_orbitals, num_particles, preserve_spin=False
        # )
        return singles + doubles

    if isinstance(integ.algo.solver, StatefulAdaptVQE):
        algo = GroundStateEigensolver(integ.algo.qubit_mapper, integ.algo.solver.solver)
    logger.info(
        "Creating QEOM"
        )

    qeom = QEOM(
        algo,
        estimator,
        excitations=my_generator,
        aux_eval_rules=EvaluationRule.ALL,
        tol=1e-6
    )

    logger.info(
        "Removing the ElectronicDensity property before starting QEOM"
    )
    problem.properties.electronic_density = None

    excited_state_result = qeom.solve(problem)

    # Print clear separation for results
    summary = f"""
{'='*80}
                         QUANTUM CALCULATION RESULTS
{'='*80}

CONFIGURATION:
-------------
Backend: {'Amazon Braket Local Simulator with ADAPT-VQE' if args.adapt else 'Amazon Braket Local Simulator'}
Number of alpha electrons: {num_alpha}
Number of beta electrons: {num_beta}
Number of orbitals: {num_orbs}
Number of shots: {estimator.options.default_shots}

CALCULATION RESULTS:
------------------
Ground State Energy: {excited_state_result.groundstate_energy if hasattr(excited_state_result, 'groundstate_energy') else 'N/A'}

Excitation Energies:
{excited_state_result.raw_result.excitation_energies}

Transition Amplitudes:
"""
    for key, values in excited_state_result.raw_result.transition_amplitudes.items():
        summary += f"\nState {key}:\n"
        for name, val in values.items():
            summary += f"  {name}: {val[0]:.6f}\n"

    summary += f"""
{'='*80}
                         CALCULATION COMPLETED
{'='*80}
Total Iterations: {solver._eval_count if hasattr(solver, '_eval_count') else 'N/A'}
Final Energy: {excited_state_result.groundstate_energy if hasattr(excited_state_result, 'groundstate_energy') else 'N/A'}
Convergence Status: {'Converged' if hasattr(solver, '_converged') and solver._converged else 'Completed'}
{'='*80}
"""

    # Print to both console and log file
    print(summary)
    logger.info(summary)

    # Save results to a JSON file
    results_dict = {
        "configuration": {
            "backend": "Amazon Braket Local Simulator",
            "method": "ADAPT-VQE" if args.adapt else "VQE",
            "num_alpha": num_alpha,
            "num_beta": num_beta,
            "num_orbs": num_orbs,
            "num_shots": estimator.options.default_shots
        },
        "results": {
            "ground_state_energy": float(excited_state_result.groundstate_energy) if hasattr(excited_state_result, 'groundstate_energy') else None,
            "excitation_energies": excited_state_result.raw_result.excitation_energies.tolist(),
            "transition_amplitudes": {
                str(key): {str(name): float(val[0]) for name, val in values.items()}
                for key, values in excited_state_result.raw_result.transition_amplitudes.items()
            }
        }
    }

    with open('quantum_calculation_results.json', 'w') as f:
        json.dump(results_dict, f, indent=4)

    logger.info("Results have been saved to 'quantum_calculation_results.json'")

class CP2KIntegration:
    def __init__(self, algo):
        self.algo = algo
        self.socket = None

    def connect_to_socket(self, host, port, unix):
        try:
            if unix:
                self.socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
                self.socket.connect(host)
            else:
                self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                self.socket.connect((host, port))
            logger.info("Connected to socket successfully.")
        except socket.error as e:
            logger.error(f"Socket connection error: {e}")
            self.socket = None

    def run(self):
        if not self.socket:
            logger.error("No socket connection available.")
            return

        try:
            # Ensure the socket is connected before sending data
            self.socket.send(self.Messages.HAVEDATA.value)
        except BrokenPipeError:
            logger.error("Broken pipe error: Unable to send data. The connection may have been closed.")
            self.reconnect()

    def reconnect(self):
        logger.info("Attempting to reconnect to the socket...")
        time.sleep(5)  # Wait before attempting to reconnect
        self.connect_to_socket(HOST, PORT, UNIX)
