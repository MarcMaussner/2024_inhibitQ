import pennylane as qml
from braket.aws import AwsDevice
from braket.devices import LocalSimulator

device_arn = "arn:aws:braket:::device/quantum-simulator/amazon/sv1"
s3_folder = ("amazon-braket-us-east-1-025066247725", "sv1-hybridjob")

dev = qml.device("braket.aws.qubit", device_arn=device_arn, s3_destination_folder=s3_folder, wires=2)

@qml.qnode(dev)
def circuit(x, y):
    qml.RX(x, wires=0)
    qml.RY(y, wires=1)
    qml.CNOT(wires=[0, 1])
    return qml.expval(qml.PauliZ(0)), qml.expval(qml.PauliZ(1))

result = circuit(0.5, 0.8)
print(f"Circuit output: {result}")

print("Circuit structure:")
print(qml.draw(circuit)(0.5, 0.8))

grad_fn = qml.grad(circuit)
grad = grad_fn(0.5, 0.8)
print(f"Gradient: {grad}")