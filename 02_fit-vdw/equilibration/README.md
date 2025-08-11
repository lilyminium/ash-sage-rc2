# Equilibrating properties

This directory contains scripts for equilibrating properties for training and validation. This step is necessary before running a refit in `../refit/`.

First, equilibration options are set-up and written out in `write-options.py`. A copy is uploaded in `options.json`.

An example of the code used to equilibrate the training data is present in `run-equilibrate-slurm.hpc3`, which calls `equilibrate-slurm.py`. A full environment is present in `conda-env.yaml`.