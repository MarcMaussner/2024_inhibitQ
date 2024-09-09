# algo.py or main.py
import pennylane as qml
from pennylane import numpy as np
from braket.jobs.metrics import log_metric
from braket.jobs import save_job_result

def qubit_rotation(num_steps=10, stepsize=0.5):
    device = qml.device("default.qubit", wires=1)

    @qml.qnode(device)
    def circuit(params):
        qml.RX(params[0], wires=0)
        qml.RY(params[1], wires=0)
        return qml.expval(qml.PauliZ(0))

    opt = qml.GradientDescentOptimizer(stepsize=stepsize)
    params = np.array([0.5, 0.75])

    for i in range(num_steps):
        params = opt.step(circuit, params)
        expval = circuit(params)
        log_metric(metric_name="expval", iteration_number=i, value=expval)

    return params.tolist()

def main():
    result = qubit_rotation()
    save_job_result({"final_params": result})
    print(f"Final parameters: {result}")

if __name__ == "__main__":
    main()