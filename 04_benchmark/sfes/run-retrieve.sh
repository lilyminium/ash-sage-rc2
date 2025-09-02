#!/bin/bash

# conda activate pontibus-alchemiscale-022

mkdir -p logs/retrieve

python retrieve.py -i scoped-keys/fsolv_openff-2.3.0rc1_key.dat -ff "Sage 2.3.0rc1" -d "FreeSolv" -o results/freesolv/openff-2.3.0rc1.csv > logs/retrieve/fsolv-openff-2.3.0rc1.log
python retrieve.py -i scoped-keys/mnsol_openff-2.3.0rc1_key.dat -ff "Sage 2.3.0rc1" -d "MNSol" -o results/mnsol/openff-2.3.0rc1.csv > logs/retrieve/mnsol-openff-2.3.0rc1.log
