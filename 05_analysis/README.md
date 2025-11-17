# Procedure

## Computing statistics for datasets (e.g. phys-prop and sfes)

1. Run `label-with-checkmol.py` and `label-with-forcefield.py` to get labels for unique SMILES in a particular dataset (see `run-label.sh`). This writes out data to `labels/`.
2. Run `run-map-labels-to-data.sh` to map those labels to each property in a dataset. What this does is take the output of benchmark datasets in `../04_benchmark` and add labels from step 1 to each property. This writes out `labelled-data.csv` in `output/*/`.
3. Run `run-compute-statistics.sh` to compute statistics for each label type.
4. Run `run-plot-aggregated-statistics.sh` to plot aggregated statistics.
5. Run `plot-labelled-scatter.py` to plot scatter plots instead of stats.
6. Run `run-get-largest-differences.sh` to get the largest meaningful differences in a dataset. Results are written to logs.
7. Also, `run-get-overrepresented-groups.sh` gets groups that are over-represented in a high-error set. This retreads some of the same code as earlier scripts and saves a `labelled-properties.csv`. 