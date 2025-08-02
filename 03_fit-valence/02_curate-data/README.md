# Data curation

The scripts in this directory select optimizations and torsiondrive records for fitting valence data.

Two steps have to occur first. 

Firstly, data must already be downloaded and preprocessed using the scripts in the `../../01_download-data/qm` directory. This is so we can easily re-filter and re-process without having to communicate with the QCArchive each time.

Note: I've uploaded my pre-processed data into that directory so you don't have to re-run it, but please make sure to **exclude the OpenFF Benchmark dataset** or else you will include test data in the training set. You can exclude the OpenFF benchmark set either by excluding the JSON file with the `-xd` flag or by deleting the `.parquet` file.

Secondly, data must be labelled using the scripts in the `label/` directory -- please see README there for more.

In my experience each selection script takes 30-60 minutes to run depending on how many parallel processes you use. Much of this is inefficient filtering that can be commented out if you're in a rush.

The expected outputs are:

- output/*.json files -- containing the training set torsions or optimizations
- counts/*.json files -- containing how many molecules matched per parameter

Please see `run-select.sh` for example of use. Each script also has a short docstring.