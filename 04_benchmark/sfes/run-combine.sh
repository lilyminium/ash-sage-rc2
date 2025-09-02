#!/bin/bash

mkdir -p logs/combine

python combine-results.py \
    -i results/freesolv   \
    -o output/freesolv.csv > logs/combine/freesolv.log

python combine-results.py \
    -i results/mnsol      \
    -o output/mnsol.csv > logs/combine/mnsol.log
