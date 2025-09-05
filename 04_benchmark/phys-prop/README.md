# Benchmarking

The code in this directory benchmarks a single property with a single force field.
While normally you would expect to benchmark an entire dataset at a time, the benchmark is set up this way to make efficient use of pre-emptible resources on a cluster. As a result, each property is saved separately in a single JSON file and needs to be collated in a final step.


The script in this directory makes use of the PreequilibratedSimulation workflow that starts from a pre-equilibrated box. As such, properties must first be equilibrated (as in `../../02_fit-vdw/equilibration`). The `stored_data` directory should either be copied to this directory, or the script should set it as an option `benchmark.py -s <path_to_stored_data>`.

Once all benchmarks are done, run `run-gather-benchmarks.sh` to save a single CSV.

## Remapping

Note -- there's some helper scripts here to re-map and copy-over pre-existing data from a different directory. That can be ignored.