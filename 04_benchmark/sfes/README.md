# SFEs with pontibus

This directory contains scripts and outputs from running and retrieving SFEs with Pontibus.

Directions:

- `run-generate.sh` runs `gen-test-[hfe/sfe]-network.py` to generate FreeSolv and MNSol networks
- `run-submit.sh` runs `submit.py` to submit said networks to Alchemiscale and write the scope keys to file.
- Run `kubectl apply -f jobs.yaml` to run workers on the Nautilus Research Platform. Monitor GPU usage here: resources-namespace-gpus?orgId=1&from=now-30m&to=now&timezone=browser&var-namespace=openforcefield&refresh=30s
- Run `run-monitor.sh` to monitor ongoing progress in terms of completed calculations.
- Run `run-retrieve.sh` to retrieve calculation results from the server.
- (If not run already): run `run-convert-data.sh` to convert experimental reference data from original sources to processed dataframes.
- Run `run-combine.sh` to combine experimental results into output files.
