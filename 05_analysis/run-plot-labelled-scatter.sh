#!/usr/bin/env bash

mkdir -p logs/plot

# # physical properties
# python plot-labelled-scatter.py                             \
#     -i  output/phys-prop/density-all/labelled-data.csv    \
#     -s  output/phys-prop/density-all/stats/Sage-2.0.0.csv                  \
#     -s  output/phys-prop/density-all/stats/Sage-2.2.1.csv                  \
#     -s  output/phys-prop/density-all/stats/Sage-2.3.0rc1.csv                  \
#     -s  output/phys-prop/density-all/stats/v1-k100.csv                  \
#     -s  output/phys-prop/density-all/stats/v3-k100.csv                  \
#     -im images/phys-prop-scatters/densities-refit                \
#     -cy     value       -cx     reference_value             \
#     -cye    uncertainty -cxe    reference_uncertainty       \
#     -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'Sage 2.3.0rc1'  -ff 'v1-k100' -ff 'v3-k100'   \
#     -u '[g/mL]' -h 4.5 -a 0.8 \
#     --plot-all-groups > logs/plot/densities-refit.log

# python plot-labelled-scatter.py                             \
#     -i  output/phys-prop/dhmix-all/labelled-data.csv    \
#     -s output/phys-prop/dhmix-all/stats/Sage-2.0.0.csv                  \
#     -s output/phys-prop/dhmix-all/stats/Sage-2.2.1.csv                  \
#     -s output/phys-prop/dhmix-all/stats/Sage-2.3.0rc1.csv                  \
#     -s output/phys-prop/dhmix-all/stats/v1-k100.csv                  \
#     -s output/phys-prop/dhmix-all/stats/v3-k100.csv                  \
#     -im images/phys-prop-scatters/dhmix-refit                \
#     -cy     value       -cx     reference_value             \
#     -cye    uncertainty -cxe    reference_uncertainty       \
#     -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'Sage 2.3.0rc1' -ff 'v1-k100' -ff 'v3-k100'   \
#     -u '[kJ/mol]' -h 4.5 -a 0.8 \
#     --plot-all-groups > logs/plot/dhmix-refit.log

# sfes

# python plot-labelled-scatter.py                                         \
#     -i  output/sfes/mnsol/labelled-data.csv                             \
#     -s  output/sfes/mnsol/stats/Sage-2.2.1-+-ELF10.csv                  \
#     -s  output/sfes/mnsol/stats/Sage-2.2.1-+-AshGC.csv                  \
#     -s  output/sfes/mnsol/stats/Sage-2.3.0rc1.csv                       \
#     -s  output/sfes/mnsol/stats/v1-k100.csv                  \
#     -s  output/sfes/mnsol/stats/v3-k100.csv                  \
#     -im images/sfe-scatters/mnsol                   \
#     -cy     'Value (kcal / mol)'                    \
#     -cx     'Reference Value (kcal / mol)'          \
#     -cye    'Uncertainty (kcal / mol)'              \
#     -cxe    'Reference Uncertainty (kcal / mol)'    \
#     -ff 'Sage 2.2.1 + ELF10' -ff 'Sage 2.2.1 + AshGC' -ff 'Sage 2.3.0rc1'  -ff 'v1-k100' -ff 'v3-k100'   \
#     -u '[kcal/mol]' -h 4.5 -a 0.8 \
#     --plot-all-groups > logs/plot/mnsol-scatter.log

python plot-labelled-scatter.py                                         \
    -i  output/sfes/mnsol/labelled-data.csv                             \
    -s  output/sfes/mnsol/stats/Sage-2.2.1-+-ELF10.csv                  \
    -s  output/sfes/mnsol/stats/Sage-2.2.1-+-AshGC.csv                  \
    -s  output/sfes/mnsol/stats/Sage-2.3.0rc1.csv                       \
    -s  output/sfes/mnsol/stats/v1-k100.csv                  \
    -s  output/sfes/mnsol/stats/v3-k100.csv                  \
    -im images/sfe-scatters/mnsol-close                  \
    -cy     'Value (kcal / mol)'                    \
    -cx     'Reference Value (kcal / mol)'          \
    -cye    'Uncertainty (kcal / mol)'              \
    -cxe    'Reference Uncertainty (kcal / mol)'    \
    -ff 'Sage 2.2.1 + ELF10' -ff 'Sage 2.2.1 + AshGC' -ff 'v3-k100'   \
    -u '[kcal/mol]' -h 4.5 -a 0.8 \
    --plot-all-groups > logs/plot/mnsol-scatter-close.log

