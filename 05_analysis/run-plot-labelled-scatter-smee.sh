#!/usr/bin/env bash

mkdir -p logs/plot

# python plot-labelled-scatter.py                             \
#     -i  output/phys-prop/density/labelled-properties.csv    \
#     -o  output/phys-prop/density/stats.csv                  \
#     -im images/phys-prop-scatters/densities-smee                \
#     -cx     value       -cy     reference_value             \
#     -cxe    uncertainty -cye    reference_uncertainty       \
#     -ff 'Sage 2.0.0' -ff 'Sage 2.2.1'  -ff 'smee-v2'   \
#     -u '[g/mL]' -h 4.5 -a 0.8 \
#     --plot-all-groups > logs/plot/densities-smee.log

python plot-labelled-scatter.py                             \
    -i  output/phys-prop/dhmix/labelled-properties.csv    \
    -o  output/phys-prop/dhmix/stats.csv                  \
    -im images/phys-prop-scatters/dhmix-smee                \
    -cx     value       -cy     reference_value             \
    -cxe    uncertainty -cye    reference_uncertainty       \
    -ff 'Sage 2.0.0' -ff 'Sage 2.2.1'  -ff 'smee-v2'   \
    -u '[kJ/mol]' -h 4.5 -a 0.8 \
    --plot-all-groups > logs/plot/dhmix-smee.log
