#!/usr/bin/env bash

mkdir logs
mkdir intermediate
mkdir output

# python clean-data.py > logs/clean-data.log 2>&1

# python blacklist-high-viscosities.py                        \
#     -i  intermediate/thermoml-cleaned.csv                   \
#     -v  viscosities/viscosity-stats.csv                     \
#     -o  intermediate/clean-filtered-viscosities.csv         \
#     -t  0.3 \
#     > logs/blacklist-high-viscosities.log 2>&1

python filter-data-initial.py                       \
    -i intermediate/clean-filtered-viscosities.csv  \
    -o output/initial-filtered.csv                  \
    -np 8 > logs/filter-data-initial.log 2>&1

python rename-property-ids.py       \
    -i output/initial-filtered.csv  \
    -o output/renamed-filtered.csv  \
    > logs/rename-property-ids.log 2>&1
