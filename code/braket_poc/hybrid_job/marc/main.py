import os
from braket.aws import AwsDevice
from braket.circuits import Circuit
from braket.devices import LocalSimulator
from braket.jobs import save_job_result
from braket.jobs.metrics import log_metric

def run_on_device(device, name):
    results = []
    bell = Circuit().h(0).cnot(0, 1)
    for count in range(5):
        task = device.run(bell, shots=100)
        counts = task.result().measurement_counts
        results.append(counts)
        print(f"  Run {count + 1} on {name}: {counts}")
        log_metric(metric_name=f"{name}_run_{count+1}", iteration_number=count, value=counts.get("00", 0) / 100)
    return results

def main():
    print("Hybrid job started!")

    try:
        # get the used device
        device_arn = os.environ["AMZN_BRAKET_DEVICE_ARN"]
        
        if "local" in device_arn.lower():
            device = LocalSimulator()
            name = "Embedded_Simulator"
        else:
            device = AwsDevice(device_arn)
            name = "SV1_Simulator" if "sv1" in device_arn.lower() else "Quantum_Device"

        results = run_on_device(device, name)
        
        total_00_counts = sum(result.get("00", 0) for result in results)
        average_00_probability = total_00_counts / len(results) / 100

        save_job_result({
            f"{name}_results": results,
            "average_00_probability": average_00_probability
        })

        print("Hybrid job completed successfully!")
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        save_job_result({"error": str(e)})

if __name__ == "__main__":
    main()