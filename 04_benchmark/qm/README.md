# QM benchmarks

This directory contains scripts and outputs for running QM benchmarks.

Firstly, data needs to be downloaded. This is demonstrated in `run-download.sh` with `download-industry-set.py` and `get-optimization/torsiondrive-data.py`.

Secondly, optimizations and torsiondrives are optimized and benchmarked in `benchmark-mm-optimization.py` and `benchmark-mm-torsion.py`, and data is written out into PyArrow datasets. From there, RMSDs and TFDs of optimized conformers from the original QM structure are calculated as in `get-rmsds-and-tfds.py`. In `get-all-to-all-rmsd.py` the lowest RMSD from every MM conformer to any QM conformer is computed, for the all-to-all ddEs calculated in `get-all-to-all-dde.py`. Lastly, all topology values for each optimized conformer are calculated with `compute-topology-values.py`. `run-benchmark.sh` executes all Python scripts for each force field.

`run-get-dde.sh` computes ddEs with `get-all-to-all-dde.py`.

`run-plot.sh` plots the images shown in `images/`.