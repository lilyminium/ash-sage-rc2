#!/usr/bin/env bash

mkdir -p logs/plot

# python plot-labelled-scatter.py                             \
#     -i  output/phys-prop/density/labelled-properties.csv    \
#     -o  output/phys-prop/density/stats.csv                  \
#     -im images/phys-prop-scatters/densities-refit                \
#     -cy     value       -cx     reference_value             \
#     -cye    uncertainty -cxe    reference_uncertainty       \
#     -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'Sage 2.3.0rc1'  -ff 'v1-k100' -ff 'v3-k100'   \
#     -u '[g/mL]' -h 4.5 -a 0.8 \
#     --plot-all-groups > logs/plot/densities-refit.log

# python plot-labelled-scatter.py                             \
#     -i  output/phys-prop/dhmix/labelled-properties.csv    \
#     -o  output/phys-prop/dhmix/stats.csv                  \
#     -im images/phys-prop-scatters/dhmix-refit                \
#     -cy     value       -cx     reference_value             \
#     -cye    uncertainty -cxe    reference_uncertainty       \
#     -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'Sage 2.3.0rc1' -ff 'v1-k100' -ff 'v3-k100'   \
#     -u '[kJ/mol]' -h 4.5 -a 0.8 \
#     --plot-all-groups > logs/plot/dhmix-refit.log

python plot-aggregated-statistics.py \
    -i output/phys-prop/density/stats.csv \
    -im images/compare-physprop/compare-refit/densities \
    -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'Sage 2.3.0rc1' -ff 'v1-k100' -ff 'v3-k100' \
    -u '[g/mL]' > logs/plot/aggregated-densities-refit.log

python plot-aggregated-statistics.py \
    -i output/phys-prop/dhmix/stats.csv \
    -im images/compare-physprop/compare-refit/dhmixes \
    -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'Sage 2.3.0rc1' -ff 'v1-k100' -ff 'v3-k100' \
    -u '[kJ/mol]' > logs/plot/aggregated-dhmixes-refit.log