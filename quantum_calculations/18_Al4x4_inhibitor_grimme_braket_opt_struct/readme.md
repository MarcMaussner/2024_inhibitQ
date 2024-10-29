# Following changes were made to integrated braket
## Changes in client-vqe-ucc.py
### Added imports:
from braket.aws import AwsDevice
from braket.devices import Devices
from braket.jobs import hybrid_job, save_job_result
from qiskit.primitives import BackendEstimator

from qiskit_braket_provider import BraketProvider
from qiskit_braket_provider import BraketLocalBackend

### Added parameter
    parser.add_argument("--braket", action="store_true")
### Use paramter to define braket estimator (for now shots = 10)
    if args.braket:
        # Configure the SV1-Braket backend with desired shots
        backend = BraketProvider().get_backend("SV1")
        estimator = BackendEstimator(backend=backend, options={"shots": 10})


## Changes in run.sh
Added "--braket" to call of client-vqe-ucc.py
