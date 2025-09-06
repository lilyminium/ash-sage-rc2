#!/bin/bash

# conda activate pontibus-alchemiscale-022

mkdir networks logs

INPUT_FILE="mnsol_freesolv_highrmse_openff-gnn-am1bcc-0.1.0-rc.3.pt_openff-nagl.sdf"

FFNAME="fb-fit-v3-single-mean-k100"
FORCEFIELD="../forcefields/${FFNAME}.offxml"


python gen-test-hfe-network.py                  \
    -i      $INPUT_FILE                         \
    -c      input-data/sage-fsolv-test-v1.csv   \
    -ff     $FORCEFIELD                         \
    -o      networks/fsolv-${FFNAME}-network.json > logs/gen-fsolv-${FFNAME}.log


python gen-test-sfe-network.py                  \
    -i      $INPUT_FILE                         \
    -c      input-data/sage-mnsol-test-v1.csv   \
    -ff     $FORCEFIELD                         \
    -o      networks/mnsol-${FFNAME}-network.json > logs/gen-mnsol-${FFNAME}.log
