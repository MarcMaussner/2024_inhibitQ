import sys

def generate_al_slab_input():
    return """&GLOBAL
  PROJECT Al_slab
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 8.0 8.0 20.0
    &END CELL
    &COORD
      Al     0.0000000000    0.0000000000    0.0000000000
      Al     2.8637069723    0.0000000000    0.0000000000
      Al     1.4318534862    2.4803039410    0.0000000000
      Al     1.4318534862    0.8267679803    2.3386838639
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
"""

def generate_al_slab_inhibitor_input():
    return """&GLOBAL
  PROJECT Al_slab_inhibitor
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 8.0 8.0 25.0
    &END CELL
    &COORD
      Al     0.0000000000    0.0000000000    0.0000000000
      Al     2.8637069723    0.0000000000    0.0000000000
      Al     1.4318534862    2.4803039410    0.0000000000
      Al     1.4318534862    0.8267679803    2.3386838639
      O      1.4318534862    0.8267679803    5.0000000000
      H      0.5000000000    0.8267679803    5.3000000000
      H      2.3637069723    0.8267679803    5.3000000000
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
"""

def generate_inhibitor_input():
    return """&GLOBAL
  PROJECT inhibitor
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
    &END CELL
    &COORD
      O      0.0000000000    0.0000000000    0.0000000000
      H      0.9318534862    0.0000000000    0.3000000000
      H     -0.9318534862    0.0000000000    0.3000000000
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
"""

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python generate_input.py <input_type>")
        print("input_type can be: al_slab, al_slab_inhibitor, or inhibitor")
        sys.exit(1)

    input_type = sys.argv[1]
    if input_type == "al_slab":
        input_content = generate_al_slab_input()
        filename = "Al_slab.inp"
    elif input_type == "al_slab_inhibitor":
        input_content = generate_al_slab_inhibitor_input()
        filename = "Al_slab_inhibitor.inp"
    elif input_type == "inhibitor":
        input_content = generate_inhibitor_input()
        filename = "inhibitor.inp"
    else:
        print("Invalid input type. Choose al_slab, al_slab_inhibitor, or inhibitor.")
        sys.exit(1)

    with open(filename, "w") as f:
        f.write(input_content)
    print(f"Generated {filename}")