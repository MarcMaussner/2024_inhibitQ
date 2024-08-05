This folder is an attempt to make the literature work in [https://archive.materialscloud.org/record/2024.66](https://archive.materialscloud.org/record/2024.66) work.


## Running Local CP2K Calculations

To run local CP2K calculations using the `run_cp2k.sh` script, follow these steps:

1. Ensure you have the necessary dependencies installed:
   - Bash
   - MPICH
   - CP2K

2. If MPICH is not installed, the script will attempt to install it using Homebrew.

3. If CP2K is not installed, the script will attempt to install it using Homebrew.

4. Run the script with the input file as an argument:
   ```sh
   ./run_cp2k.sh <input_file>
   ```

Replace <input_file> with the path to your CP2K input file.

The script will execute CP2K with the provided input file and save the output to a file with the same base name but with a .out extension.

If the run is successful, the output and restart files will be moved to a directory named after the input file (without the extension).

Example:
```sh
./run_cp2k.sh H2O-dft-ls.NREP2.inp
```
This will run CP2K with the H2O-dft-ls.NREP2.inp input file and save the output to H2O-dft-ls.NREP2.out.

If you are using ML potentials, you should follow these steps:

To use these scripts:

For MatGL:
```sh
pip install matgl
install torch torchvision
pip install dgl -f https://data.dgl.ai/wheels/repo.html
pip install git+https://github.com/materialsvirtuallab/matgl.git
```
For CHGNet:
```sh
pip install chgnet
pip install git+https://github.com/CederGroupHub/chgnet
pip install dftd3
```
Run the scripts:
```sh
python matgl_script.py
python chgnet_script.py
```

with GPAW
```sh
brew install libxc
export CFLAGS="-I$(brew --prefix libxc)/include"
export LDFLAGS="-L$(brew --prefix libxc)/lib"
pip install gpaw
 gpaw install-data
 ```

 To redo the examples in the literature publication, you can check the literature folder then MgO