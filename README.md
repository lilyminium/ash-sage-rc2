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