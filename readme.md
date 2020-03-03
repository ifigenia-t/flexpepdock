# Download and install Rosetta

```
# Obtain an academic license for Rosetta and PyRosetta

pip3 search PyRosetta
unzip rosetta_bin_mac_3.11_bundle.tgz
tar xvf rosetta_bin_mac_3.11_bundle.tgz
tar xvf PyRosetta4.Release.python37.mac.release-242.tar.bz2
cd rosetta_bin_mac_2019.35.60890_bundle
cd ..
cd PyRosetta4.Release.python37.mac.release-242
cd main/source

# Rosetta compile
./scons.py -j 2 mode=release bin
cd ..
cd Downloads/flex

cd setup && sudo python setup.py install
import pyrosetta;pyrosetta.int()
python 3

```


# FlexPepDock


```
python3 -m venv venv
source venv/bin/activate
cd PyRosetta4.Release.python37.mac.release-242/setup
pip install .
cd flex
python3 ./run_fpb.py -p 3twr_clean.pdb -l test_peptide.txt -c 8

```

brew install slurm


```
# Clean the Score files before each run
./clean.sh

# Run the programme
./run_fpb_new.py -p 3twr_clean.pdb -l test_peptide.txt -cid H -rcid D
```