## Synopsis

**toolbx_pdb** is a set of tools for manipulation and execution of tasks on
PDB file ensembles.

## Examples

Below are example of uses cases for the toolbx_pdb scripts.

### Binding pocket conformation ensemble analysis

**Generate an interaction fingerprints (IFP) diagram**
Generate an IFP diagram for all binding pockets in "complex_ensemble"
directory. The complexes are ligand bound proteins. Each residue that interact
with its bound ligand is listed, and the combined list of all residues is used
to generate the IFP. A custom set of interaction types can be defined with
the -customFprint flag.
```
pdb_confEns.py 'IFP_example' 'complex_ensemble' -customFprint '11111111000'
```

**Cluster complexes by IFP**
Cluster the complexes based on their IFP and generate a dendrogram of that
clustering. An additional template complex is added for comparison.
```
pdb_confEns.py 'IFP_cluster_example' 'complex_ensemble' -templatePath 'structures/pocket1.pdb' -dendro
```

**Analyze binding pocket ensemble conformation with principal component analysis (PCA)**
Analyze the range of motion across the binding pocket ensemble with PCA on
all binding pocket carbon alphas. A template and additional binding pocket
conformation are added here for comparisons.
```
pdb_confEns.py 'PCA_example' 'complex_ensemble' -templatePath 'structures/pocket1.pdb' -additionalPaths 'structures/pocket2.pdb' -pca
```

### Batch protein superimposition
Superimpose all proteins from the directory "complex_ensemble" onto the
structure of "structures/pocket1.pdb"
```
pdb_superimpose.py "complex_ensemble" "structures/pocket1.pdb"
```

### Batch add ligand charge
Guess the charge on the bound ligand by updating its protonation state. Proper
treatment of the ligand charge is necessary for IFP calculations.
```
pdb_addLigand.py "complex_ensemble"
```

### Batch get B-Factor
Display B-factor data from X-ray structure PDB files.
```
pdb_bfac.py "pdb_structures"
```

## Motivation

This set of tools simplifies the analysis of a ligand-bound binding pocket
ensemble. IFPs capture information of both the ligand and binding pocket
conformation by focusing on their interaction pattern. PCA on binding pocket
residue captures the main collective movements in a simple visual output. These
scripts were developed by Thomas Coudrat during his PhD.

## Contributions
Feel free to create an issue or submit a pull request. You can also contact me
(Thomas Coudrat) if you have questions related to this project.

## Installation
* Install Anaconda for Python 3.5 (required)
* Install OpenEye OEChem toolkit (OpenEye-python3-ubuntu-12.04-x64-2014.10.2) (optional: for interaction fingerprints)
    * set path: export PYTHONPATH="path/to/OpenEye-python3-ubuntu-12.04-x64-2014.10.2"
    * Requires software license
* Install ICM 3.8-4 (optional: for superimposition)
    * set environment variable: export ICMHOME="path/to/icmDirectory"
    * Requires software license

## License
This project is licensed under the MIT license
