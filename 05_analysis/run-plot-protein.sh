#!/usr/bin/env bash

mkdir -p logs/plot

# SFE_FF_KWARGS=(
#     -ff 'Sage 2.2.1 + ELF10'
#     -ff 'Sage 2.2.1 + AshGC'
#     -ff 'Sage 2.3.0rc1'
#     -ff 'v1-k100'
#     -ff 'v3-k100'
# )

# python plot-aggregated-statistics.py    \
#     "${SFE_FF_KWARGS[@]}"               \
#     -i output/sfes/freesolv/stats/Sage-2.2.1-+-AmberTools.csv   \
#     -i output/sfes/freesolv/stats/Sage-2.2.1-+-ELF10.csv   \
#     -i output/sfes/freesolv/stats/Sage-2.2.1-+-AshGC.csv   \
#     -i output/sfes/freesolv/stats/Sage-2.3.0rc1.csv   \
#     -i output/sfes/freesolv/stats/v1-k100.csv   \
#     -i output/sfes/freesolv/stats/v3-k100.csv   \
#     -im images/compare-refit/freesolv-stats \
#     -u '[kcal/mol]' > logs/plot/aggregated-stats_freesolv-refit.log

# python plot-aggregated-statistics.py    \
#     "${SFE_FF_KWARGS[@]}"               \
#     -i output/sfes/mnsol/stats/Sage-2.2.1-+-AmberTools.csv   \
#     -i output/sfes/mnsol/stats/Sage-2.2.1-+-ELF10.csv   \
#     -i output/sfes/mnsol/stats/Sage-2.2.1-+-AshGC.csv   \
#     -i output/sfes/mnsol/stats/Sage-2.3.0rc1.csv   \
#     -i output/sfes/mnsol/stats/v1-k100.csv   \
#     -i output/sfes/mnsol/stats/v3-k100.csv   \
#     -im images/compare-refit/mnsol-stats \
#     -u '[kcal/mol]' > logs/plot/aggregated-stats_mnsols-refit.log

# python plot-aggregated-statistics.py    \
#     "${SFE_FF_KWARGS[@]}"               \
#     -i output/sfes/tfes/stats/Sage-2.2.1-+-AmberTools.csv   \
#     -i output/sfes/tfes/stats/Sage-2.2.1-+-ELF10.csv   \
#     -i output/sfes/tfes/stats/Sage-2.2.1-+-AshGC.csv   \
#     -i output/sfes/tfes/stats/Sage-2.3.0rc1.csv   \
#     -i output/sfes/tfes/stats/v1-k100.csv   \
#     -i output/sfes/tfes/stats/v3-k100.csv   \
#     -im images/compare-refit/tfes-stats \
#     -u '[kcal/mol]' > logs/plot/aggregated-stats_tfes-refit.log


# # == physprop kwargs ==

# PHYS_PROP_FF_KWARGS=(
#     -ff 'Sage 2.0.0'
#     -ff 'Sage 2.2.1'
#     -ff 'Sage 2.3.0rc1'
#     -ff 'v1-k100'
#     -ff 'v3-k100'
# )

for dataset in density-all ; do
python plot-aggregated-statistics.py \
    "${PHYS_PROP_FF_KWARGS[@]}" \
    -i output/phys-prop/${dataset}/stats/Sage-2.1.0.csv \
    -i output/phys-prop/${dataset}/stats/Null-AAQAA3-3-+-OPC3.csv \
    -i output/phys-prop/${dataset}/stats/Null-4-mer-AAQAA3-2-+-OPC3.csv \
    -i output/phys-prop/${dataset}/stats/v3-k100.csv \
    -im images/compare-protein/${dataset} \
    -ff 'Sage 2.1.0' -ff 'Null 4-mer AAQAA3-2 + OPC3' -ff 'Null AAQAA3-3 + OPC3' -ff 'v3-k100' \
    -u '[g/mL]' > logs/plot/aggregated-stats_protein_${dataset}-refit.log
done
for dataset in dhmix-all ; do
python plot-aggregated-statistics.py \
    "${PHYS_PROP_FF_KWARGS[@]}" \
    -i output/phys-prop/${dataset}/stats/Sage-2.1.0.csv \
    -i output/phys-prop/${dataset}/stats/Null-AAQAA3-3-+-OPC3.csv \
    -i output/phys-prop/${dataset}/stats/Null-4-mer-AAQAA3-2-+-OPC3.csv \
    -i output/phys-prop/${dataset}/stats/v3-k100.csv \
    -im images/compare-protein/${dataset} \
    -ff 'Sage 2.1.0' -ff 'Null 4-mer AAQAA3-2 + OPC3' -ff 'Null AAQAA3-3 + OPC3' -ff 'v3-k100' \
    -u '[kJ/mol]' > logs/plot/aggregated-stats_protein_${dataset}-refit.log
done
