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


# Running the FlexPepDock Protocol


# Input Parameters

```
The program supports the following modes:
-   	Input Parameters:
o   Local .pdb file.
o   PDB ID, with which the relevant structure can be downloaded from the Protein Data Bank (https://www.rcsb.org/).
o   Number of the peptide and receptor chains.
o   Sequence of the peptide.
-   	Peptide Lists Parameters:
o   List of different peptides.
o   Mode of mutation we want to test so that respective peptide lists can be made.
§  Alanine Scanning: Mutation of every residue in the peptide in each position, to an alanine amino acid.
§  All by All Scanning: Mutation of every residue of the peptide in each position by all the other standard amino acids.
-   	Running Parameters:
o   Minimization including receptor backbone argument that improves the model’s accuracy.
o   Refinement that allows the refinement of coarse peptide-protein models to near-native accuracy.
o   Both Minimization and Refinement parameters

````

# Workflow 

```
These are the basic steps of the FlexPepDock Workflow:
PDB file is downloaded or located and a list of mutated peptides is created from the original peptide depending on mutation mode.
For each peptide in the list:
The peptide is threaded onto the template peptide backbone by using the fixbb design from the Rosetta algorithms (Leaver-Fay et al., 2011).
The threaded structure is optimized by using the FlexPepDock protocols by running either the minimization or/and the refinement.
We repeat the refinement and/or minimization protocols, according to the arguments, 5 or 100 times.
Create the score file.
For every model the reweighted score is calculated as following:
Reweighted score = Total score + Peptide score + Interface score
Where:
Total score = the total score of the complex
Peptide score = the internal energy of the peptide
Interface score = the sum of the energy contributions of the interface residues from the peptide and the receptor
 
 
Normalise the PSSM by: (pep_score - ala_score) / (base_score - ala_score)
                                	where: pep_score the positional score for each peptide
                                            	ala_score: the score of a peptide made entirely of alanines with length as the original peptide
                                            	base_score: the score of the original peptide
5. 	Export the PSSM.


```

```
brew install slurm

# Clean the Score files before each run
./clean.sh

# Run the programme
./run_fpb_new.py -p 3twr_clean.pdb -l test_peptide.txt -cid H -rcid D
```
