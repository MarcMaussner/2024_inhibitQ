To run the Al111 use command:

Latest working commands:
```bash
python -u client-vqe-ucc.py --nalpha 1 --nbeta 1 --norbs 5 2>&1 | tee python_output.log
docker run -v $PWD:/mnt --shm-size=1g -u $(id -u):$(id -g) -it --rm cp2k/cp2k mpiexec -genv OMP_NUM_THREADS=2 -np 18 cp2k Al111_active_space.inp
```


envorinment:
```bash
# install pyenv # skip this if you already have pyenv installed or use your own python environment
sudo apt upgrade
sudo apt install make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python3-openssl git
curl https://pyenv.run | bash
#coupy and paste the following lines to your .bashrc
export PATH="$HOME/.pyenv/bin:$PATH"
eval "$(pyenv init --path)"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
#run this command in terminal
source ~/.bashrc
pyenv install --list | grep " 3.12"
 pyenv install 3.12.5
pyenv global 3.12.5
pip install qiskit qiskit-aer qiskit-nature
git clone -b fix-qiskit1.x-compatibility https://github.com/KarimElgammal/qiskit-nature-cp2k.git
cd qiskit-nature-cp2k
pip install -e .
```

```bash
pip list
Package            Version     Editable project location
------------------ ----------- -------------------------------
dill               0.3.8
h5py               3.11.0
mpmath             1.3.0
numpy              2.1.0
pbr                6.0.0
pip                24.2
psutil             6.0.0
python-dateutil    2.9.0.post0
qiskit             1.2.0
qiskit-aer         0.14.2
qiskit-algorithms  0.3.0
qiskit-nature      0.7.2
qiskit_nature_cp2k 0.0.1       /home/ubuntu/qiskit-nature-cp2k
rustworkx          0.15.1
scipy              1.14.1
setuptools         73.0.1
six                1.16.0
stevedore          5.2.0
symengine          0.11.0
sympy              1.13.2
typing_extensions  4.12.2
```