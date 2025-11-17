#!/bin/bash

# conda activate pontibus-alchemiscale-022

mkdir -p logs/retrieve

# python retrieve.py -i scoped-keys/fsolv_openff-2.3.0rc1_key.dat -ff "Sage 2.3.0rc1" -d "FreeSolv" -o results/freesolv/openff-2.3.0rc1.csv > logs/retrieve/fsolv-openff-2.3.0rc1.log
# python retrieve.py -i scoped-keys/mnsol_openff-2.3.0rc1_key.dat -ff "Sage 2.3.0rc1" -d "MNSol" -o results/mnsol/openff-2.3.0rc1.csv > logs/retrieve/mnsol-openff-2.3.0rc1.log
# python retrieve.py -i scoped-keys/fsolv_fb-fit-v1-single-mean-k100_key.dat -ff "v1-k100" -d "FreeSolv" -o results/freesolv/v1-k100.csv > logs/retrieve/fsolv-v1-k100.log
# python retrieve.py -i scoped-keys/mnsol_fb-fit-v1-single-mean-k100_key.dat -ff "v1-k100" -d "MNSol" -o results/mnsol/v1-k100.csv > logs/retrieve/mnsol-v1-k100.log
# python retrieve.py -i scoped-keys/fsolv_fb-fit-v3-single-mean-k100_key.dat -ff "v3-k100" -d "FreeSolv" -o results/freesolv/v3-k100.csv > logs/retrieve/fsolv-v3-k100.log
# python retrieve.py -i scoped-keys/mnsol_fb-fit-v3-single-mean-k100_key.dat -ff "v3-k100" -d "MNSol" -o results/mnsol/v3-k100.csv > logs/retrieve/mnsol-v3-k100.log
# python retrieve.py -i scoped-keys/fsolv_fb-fit-v0-single-mean-k20_key.dat -ff "v0-k20" -d "FreeSolv" -o results/freesolv/v0-k20.csv > logs/retrieve/fsolv-v0-k20.log
# python retrieve.py -i scoped-keys/mnsol_fb-fit-v0-single-mean-k20_key.dat -ff "v0-k20" -d "MNSol" -o results/mnsol/v0-k20.csv > logs/retrieve/mnsol-v0-k20.log
python retrieve.py -i scoped-keys/mnsol_fit-iter-1_key.dat -ff "fit-iter-1" -d "MNSol" -o results/mnsol/fit-iter-1.csv > logs/retrieve/mnsol-fit-iter-1.log
python retrieve.py -i scoped-keys/mnsol_vdw-refit_key.dat -ff "vdw-refit" -d "MNSol" -o results/mnsol/vdw-refit.csv > logs/retrieve/mnsol-vdw-refit.log
# python retrieve.py -i scoped-keys/fsolv_vdw-refit_key.dat -ff "vdw-refit" -d "FreeSolv" -o results/freesolv/vdw-refit.csv > logs/retrieve/fsolv-vdw-refit.log