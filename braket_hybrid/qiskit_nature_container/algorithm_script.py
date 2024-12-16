import os
import time
from datetime import datetime

import numpy as np

from qiskit.circuit.library import TwoLocal
from qiskit.primitives import BackendEstimator
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.minimum_eigensolvers import VQE
from qiskit_algorithms.optimizers import SLSQP

from braket.jobs import save_job_result
from braket.jobs.metrics import log_metric
from braket.tracking import Tracker

from qiskit_braket_provider import BraketProvider

def main():
    cost_tracker = Tracker().start()
    np.random.seed(42)

    """Decorated function that will be run in the docker container."""
    backend = BraketProvider().get_backend("SV1")

    h2_op = SparsePauliOp(
        ["II", "IZ", "ZI", "ZZ", "XX"],
        coeffs=[
            -1.052373245772859,
            0.39793742484318045,
            -0.39793742484318045,
            -0.01128010425623538,
            0.18093119978423156,
        ],
    )

    estimator = BackendEstimator(backend=backend, options={"shots": 10})
    ansatz = TwoLocal(rotation_blocks="ry", entanglement_blocks="cz")
    slsqp = SLSQP(maxiter=1)

    vqe = VQE(estimator=estimator, ansatz=ansatz, optimizer=slsqp)

    vqe_result = vqe.compute_minimum_eigenvalue(h2_op)

    # Save the results of the VQE computation.
    save_job_result(
        {
            "VQE": {
                "eigenvalue": vqe_result.eigenvalue.real,
                "optimal_parameters": list(vqe_result.optimal_parameters.values()),
                "optimal_point": vqe_result.optimal_point.tolist(),
                "optimal_value": vqe_result.optimal_value.real,
            }
        }
    )