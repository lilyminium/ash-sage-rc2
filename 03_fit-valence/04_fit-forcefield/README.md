# Fitting the force field

The scripts here generate the actual ForceBalance input files, and run the fit.

We use `run-setup-k100.sh` to write ForceBalance input files with looser priors than were used in the Sage 2.2.1 fit. I recommend you read through the script as there are some workarounds that are fairly specific to our software stack and computing environment. For example, it uses a different Python environment to generate the input files than to actually run the fit. It also temporarily edits force field files to get around a software bug, and automatically copies files across different file systems.

The script `runfiles/run-fit.sh` actually execute the fit. Again, I recommend reading through this script. The main important part is calling `ForceBalance.py optimize.in`; the rest is a lot of specific tricks for our computing environment (e.g. copying across the conda environment, working around permissions and quota issues, and so on).
