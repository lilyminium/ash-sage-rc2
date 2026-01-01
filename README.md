# ash-sage-rc2

This is a **working directory** for the Sage 2.3.0 fit.

Sage 2.3.0 will be the first OpenFF force field with a neural network charge model.
The scripts here carry out both a vdW and a valence fit to the charge model.

Steps are broken up into the following directories:

- 01_download-data: downloading training and benchmarking data of various types (e.g. physical properties, QM). 
- 02_fit-vdw: Fitting vdW parameters using physical property data and OpenFF Evaluator.
- 03_fit-valence: Fitting valence parameters using QM data and ForceBalance.
- 04_benchmark: handling QM, physical property, SFE, RBFE benchmarks.
- 05_analysis: analysis and plotting
- 06_additional-analysis: some further investigation into sfes
- 07_plot-for-preprint: some plotting code for generating images for a preprint. This code is liable to change.


In general we provide both Python scripts for executing steps and shell scripts demonstrating execution in each directory. Where possible files and data are provided, but in some cases size limits or licensing does not permit this.

## Usage

This force field will get released as Sage 2.3.0. When it is out, you will be able to use the force field as with any other OpenFF force field:

```python
from openff.toolkit import Molecule, ForceField
ff = ForceField("openff-2.3.0.offxml")
caffeine = Molecule.from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
interchange = ff.create_interchange(caffeine.to_topology())
openmm_system = interchange.to_openmm()
```

Alternatively, you can use AshGC ('openff-gnn-am1bcc-1.0.0') directly to assign charges:

```python
from openff.toolkit import Molecule
caffeine = Molecule.from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
caffeine.assign_partial_charges("openff-gnn-am1bcc-1.0.0.pt")
```
