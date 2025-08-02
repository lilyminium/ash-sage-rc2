#!/bin/bash

conda activate evaluator-050

DIRECTORY="../../../01_download-data/qm/data/tables"

python filter-seesaw-sx4.py         \
    -dt optimization                \
    -np 8                           \
    -i $DIRECTORY/optimization      \
    -o failed-seesaw-sx4.dat

python filter-seesaw-sx4.py         \
    -dt torsiondrive                \
    -np 8                           \
    -i $DIRECTORY/torsiondrive      \
    -o failed-seesaw-sx4-tds.dat
