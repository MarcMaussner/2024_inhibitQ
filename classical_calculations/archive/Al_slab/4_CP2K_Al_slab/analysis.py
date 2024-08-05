import os

def analyze_output(filename):
    energy = None
    with open(filename, 'r') as f:
        for line in f:
            if "ENERGY| Total FORCE_EVAL" in line:
                energy = float(line.split()[-1])
                break
    return energy

def main():
    outputs = ['Al_slab-1.out', 'Al_slab_inhibitor-1.out', 'inhibitor-1.out']
    energies = {}

    for output in outputs:
        if os.path.exists(output):
            energy = analyze_output(output)
            energies[output] = energy
            print(f"Energy for {output}: {energy} a.u.")
        else:
            print(f"Output file {output} not found.")

    if len(energies) == 3:
        binding_energy = energies['Al_slab_inhibitor-1.out'] - (energies['Al_slab-1.out'] + energies['inhibitor-1.out'])
        print(f"Binding Energy: {binding_energy} a.u.")

if __name__ == "__main__":
    main()