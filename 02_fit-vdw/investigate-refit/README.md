# Investigate refit

The scripts here analyse the refit after completion. See run.sh for what was run.

Collate the training data (ensure you have downloaded `optimize.tmp/` into `../refit`.) The script assumes you have run for 15 iterations of training; you may need to modify this as a hardcoded parameter if that is not the case.

```
python collate-data.py -i ../refit -o output/training-per-iteration.csv
```

Quickly compare RMSEs of properties from the initial FF to the last iteration with `compare-iterations-from-training[-N].py`. The `-N` accepts a substructure pattern parameter to only show properties that match a SMARTS pattern.

**Note: this is a quick short-hand comparison. The last iteration is not necessarily the optimized force field. For true comparison of training data, see `../../04_benchmark/physprop` for scripts to benchmark properly with multiple replicates per property.**

Look at parameter progress over iterations with `compare-parameters-over-iterations.py`. This is interesting to look at whether parameters are generally converged around a particular average, or still trending upwards or downwards.

Finally, use `plot-parameter-changes-percent.py` to plot the % change per parameter. Note, some parameters (e.g. figure size, x-axis limits) are hardcoded and may need to be adjusted to fit your data.
