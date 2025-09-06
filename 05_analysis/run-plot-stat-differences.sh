#!/usr/bin/env bash

# python plot-stat-differences.py \
#     -i output/phys-prop/density/stats.csv \
#     -o images/compare-physprop/compare-smee/densities

# python plot-stat-differences.py \
#     -i output/phys-prop/dhmix/stats.csv \
#     -o images/compare-physprop/compare-smee/dhmixes

python plot-stat-differences.py \
    -t "v1-k100" \
    -i output/phys-prop/density/stats.csv \
    -o images/compare-physprop/compare-v1/densities

python plot-stat-differences.py \
    -t "v1-k100" \
    -i output/phys-prop/dhmix/stats.csv \
    -o images/compare-physprop/compare-v1/dhmixes

python plot-stat-differences.py \
    -t "v3-k100" \
    -i output/phys-prop/density/stats.csv \
    -o images/compare-physprop/compare-v3/densities

python plot-stat-differences.py \
    -t "v3-k100" \
    -i output/phys-prop/dhmix/stats.csv \
    -o images/compare-physprop/compare-v3/dhmixes

