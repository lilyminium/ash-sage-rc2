# Downloading data

The scripts here download initial data from the MolSSI QCArchive:

* optimization and torsiondrive data (for training and testing)
* Hessian data (for MSM values)

Please see `run-download.sh` for input/output options, and *.log files for logging output.


## Optimizations and Torsiondrives

The following datasets were pooled initially before further data curation. The following simple filters were applied to remove unsuitable data:

* if a record was not complete (had "COMPLETE" status)
* if connectivity changed
* if stereochemistry was unperceivable
* if there were H-bonds (for torsiondrives)
* if there were RDKit connection errors when parsing SMILES

**Optimization data**

    * OpenFF Optimization Set 1
    * SMIRNOFF Coverage Set 1
    * OpenFF VEHICLe Set 1
    * OpenFF Discrepancy Benchmark 1
    * OpenFF Ehrman Informative Optimization v0.2
    * Pfizer discrepancy optimization dataset 1
    * FDA optimization dataset 1
    * Kinase Inhibitors: WBO Distributions
    * OpenFF Gen 2 Opt Set 1 Roche
    * OpenFF Gen 2 Opt Set 2 Coverage
    * OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy
    * OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy
    * OpenFF Gen 2 Opt Set 5 Bayer
    * OpenFF Sandbox CHO PhAlkEthOH v1.0
    * OpenFF Aniline Para Opt v1.0
    * OpenFF Gen2 Optimization Dataset Protomers v1.0
    * OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0
    * OpenFF Iodine Chemistry Optimization Dataset v1.0
    * XtalPi Shared Fragments OptimizationDataset v1.0
    * XtalPi 20-percent Fragments OptimizationDataset v1.0
    * OpenFF Torsion Benchmark Supplement v1.0
    * OpenFF Torsion Multiplicity Optimization Training Coverage Supplement v1.0
    * OpenFF Torsion Multiplicity Optimization Benchmarking Coverage Supplement v1.0
    * OpenFF Iodine Fragment Opt v1.0
    * OpenFF Sulfur Optimization Training Coverage Supplement v1.0
    * OpenFF Sulfur Optimization Benchmarking Coverage Supplement v1.0
    * OpenFF Lipid Optimization Training Supplement v1.0
    * OpenFF Lipid Optimization Benchmark Supplement v1.0
    * OpenFF Cresset Additional Coverage Optimizations v4.0
    * OpenFF Protein PDB 4-mers v4.0
    * OpenFF Additional Generated ChEMBL Optimizations v4.0
  

**Torsiondrives**

    * OpenFF Group1 Torsions
    * OpenFF Group1 Torsions 2
    * OpenFF Group1 Torsions 3
    * SMIRNOFF Coverage Torsion Set 1
    * OpenFF Substituted Phenyl Set 1
    * Pfizer discrepancy torsion dataset 1
    * OpenFF Primary Benchmark 1 Torsion Set
    * OpenFF Gen 2 Torsion Set 1 Roche 2
    * OpenFF Gen 2 Torsion Set 2 Coverage 2
    * OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2
    * OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2
    * OpenFF Gen 2 Torsion Set 5 Bayer 2
    * OpenFF Gen 2 Torsion Set 6 supplemental 2
    * OpenFF Fragmenter Validation 1.0
    * OpenFF DANCE 1 eMolecules t142 v1.0
    * OpenFF Rowley Biaryl v1.0
    * OpenFF-benchmark-ligand-fragments-v1.0
    * OpenFF Protein Fragments TorsionDrives v1.0
    * OpenFF WBO Conjugated Series v1.0
    * OpenFF Amide Torsion Set v1.0
    * OpenFF-benchmark-ligand-fragments-v2.0
    * OpenFF multiplicity correction torsion drive data v1.1
    * OpenFF Protein Capped 3-mer Omega v1.0
    * XtalPi Shared Fragments TorsiondriveDataset v1.0
    * OpenFF Torsion Coverage Supplement v1.0
    * OpenFF RNA Dinucleoside Monophosphate TorsionDrives v1.0
    * XtalPi 20-percent Fragments TorsiondriveDataset v1.0
    * OpenFF Torsion Drive Supplement v1.0
    * OpenFF Torsion Multiplicity Torsion Drive Coverage Supplement v1.0
    * OpenFF Phosphate Torsion Drives v1.0
    * OpenFF Alkane Torsion Drives v1.0
    * OpenFF Cresset Additional Coverage TorsionDrives v4.0
    * OpenFF Additional Generated ChEMBL TorsionDrives 4.0
    * OpenFF Additional Generated Guanidine Amidine Derivative and O-Linker TorsionDrives 4.0
    * OpenFF Gen3 Torsion Set v1.0


Molecules containing iodine in the following datasets were removed, due to known issues:

    * OpenFF Discrepancy Benchmark 1
    * OpenFF Gen 2 Opt Set 2 Coverage
    * OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy
    * SMIRNOFF Coverage Set 1
    * OpenFF Ehrman Informative Optimization v0.2
    * FDA optimization dataset 1
    * Kinase Inhibitors: WBO Distributions
    * OpenFF Gen 2 Torsion Set 2 Coverage 2
    * OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2


## Hessians

Hessian data of all optimization data downloaded in the previous step, was also downloaded. For each molecule, the lowest energy conformer was selected and included in the dataset.