#!/bin/bash


python download-qcdata-tables.py -o data/tables > download.log 2>&1

python download-hessian-data.py                         \
    -i      data/tables/optimization                    \
    -o      hessian-data                                    \
    -x      ../../03_fit-valence/02_curate-data/bad-data/bad-qcarchive_ids.dat > download-hessian.log 2>&1
