from braket.circuits import Circuit 
from braket.devices import LocalSimulator 
device = LocalSimulator() 
bell = Circuit().h(0).cnot(0, 1) 
print(device.run(bell, shots=100).result().measurement_counts)

