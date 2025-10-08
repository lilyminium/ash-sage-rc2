#!/usr/bin/env bash

mkdir -p logs/plot

SFE_FF_KWARGS=(
    -ff 'Sage 2.2.1 + ELF10'
    -ff 'Sage 2.2.1 + AshGC'
    -ff 'Sage 2.3.0rc1'
    -ff 'v3-k100'
)

python plot-aggregated-statistics.py    \
    "${SFE_FF_KWARGS[@]}"               \
    -id output/sfes/mnsol/stats \
    -im images/compare-refit/mnsol-stats \
    -u '[kcal/mol]' > logs/plot/aggregated-stats_mnsols-refit.log
    
python plot-aggregated-statistics.py    \
    "${SFE_FF_KWARGS[@]}"               \
    -id output/sfes/freesolv/stats \
    -im images/compare-refit/freesolv-stats \
    -u '[kcal/mol]' > logs/plot/aggregated-stats_freesolv-refit.log

python plot-aggregated-statistics.py    \
    "${SFE_FF_KWARGS[@]}"               \
    -id output/sfes/tfes/stats \
    -im images/compare-refit/tfes-stats \
    -u '[kcal/mol]' > logs/plot/aggregated-stats_tfes-refit.log



# SFE_FF_KWARGS=(
#     -ff 'Sage 2.2.1 + ELF10'
#     -ff 'Sage 2.2.1 + AshGC'
#     -ff 'fit-iter-1'
#     -ff 'vdw-refit'
#     -ff 'v3-k100'
# )

# python plot-aggregated-statistics.py    \
#     "${SFE_FF_KWARGS[@]}"               \
#     -id output/sfes/mnsol/stats \
#     -im images/compare-refit-mnsol/mnsol-stats \
#     -u '[kcal/mol]' #> logs/plot/aggregated-stats_mnsols-refit.log
    
# python plot-aggregated-statistics.py    \
#     "${SFE_FF_KWARGS[@]}"               \
#     -id output/sfes/freesolv/stats \
#     -im images/compare-refit-mnsol/freesolv-stats \
#     -u '[kcal/mol]' #> logs/plot/aggregated-stats_freesolv-refit.log

# python plot-aggregated-statistics.py    \
#     "${SFE_FF_KWARGS[@]}"               \
#     -id output/sfes/tfes/stats \
#     -im images/compare-refit-mnsol/tfes-stats \
#     -u '[kcal/mol]' #> logs/plot/aggregated-stats_tfes-refit.log


# # == physprop kwargs ==

PHYS_PROP_FF_KWARGS=(
    # -ff 'Sage 2.0.0'
    # -ff 'Sage 2.2.1'
    # -ff 'Sage 2.3.0rc1'
    # -ff 'v3-k100'
)

for dataset in density-all density-validation ; do
python plot-aggregated-statistics.py \
    "${PHYS_PROP_FF_KWARGS[@]}" \
    -i output/phys-prop/${dataset}/stats/Sage-2.0.0.csv \
    -i output/phys-prop/${dataset}/stats/Sage-2.2.1.csv \
    -i output/phys-prop/${dataset}/stats/Sage-2.2.1-+-AshGC.csv \
    -i output/phys-prop/${dataset}/stats/Sage-2.3.0rc1.csv \
    -i output/phys-prop/${dataset}/stats/v0-k20.csv \
    -i output/phys-prop/${dataset}/stats/v1-k100.csv \
    -i output/phys-prop/${dataset}/stats/v3-k100.csv \
    -im images/compare-refit/${dataset} \
    -ff 'Sage 2.2.1' -ff 'Sage 2.2.1 + AshGC' -ff 'Sage 2.3.0rc1' -ff 'v3-k100' \
    -u '[g/mL]' > logs/plot/aggregated-stats_${dataset}-refit.log
done
for dataset in dhmix-all dhmix-validation ; do
python plot-aggregated-statistics.py \
    "${PHYS_PROP_FF_KWARGS[@]}" \
    -i output/phys-prop/${dataset}/stats/Sage-2.0.0.csv \
    -i output/phys-prop/${dataset}/stats/Sage-2.2.1.csv \
    -i output/phys-prop/${dataset}/stats/Sage-2.2.1-+-AshGC.csv \
    -i output/phys-prop/${dataset}/stats/Sage-2.3.0rc1.csv \
    -i output/phys-prop/${dataset}/stats/v0-k20.csv \
    -i output/phys-prop/${dataset}/stats/v1-k100.csv \
    -i output/phys-prop/${dataset}/stats/v3-k100.csv \
    -im images/compare-refit/${dataset} \
    -ff 'Sage 2.2.1' -ff 'Sage 2.2.1 + AshGC' -ff 'Sage 2.3.0rc1' -ff 'v3-k100' \
    -u '[kJ/mol]' > logs/plot/aggregated-stats_${dataset}-refit.log
done

# QM_FF_KWARGS=(
#     -ff 'Sage 2.0.0'
#     -ff 'Sage 2.2.1'
#     -ff 'Sage 2.3.0rc1'
#     -ff 'v1-k100'
#     -ff 'v3-k100'
# )


# python plot-aggregated-statistics.py \
#     "${QM_FF_KWARGS[@]}" \
#     -i output/qm/mm-qm-diff/stats/Bonds.csv \
#     -im images/compare-refit/qm-bonds \
#     -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'Sage 2.3.0rc1' -ff 'v1-k100' -ff 'v3-k100' \
#     -u '[$\AA$]' > logs/plot/aggregated-stats_bonds-refit.log

# python plot-aggregated-statistics.py \
#     "${QM_FF_KWARGS[@]}" \
#     -i output/qm/mm-qm-diff/stats/Angles.csv \
#     -im images/compare-refit/qm-angles \
#     -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'Sage 2.3.0rc1' -ff 'v1-k100' -ff 'v3-k100' \
#     -u '[°]' > logs/plot/aggregated-stats_angles-refit.log

RBFE_FF_KWARGS=(
    -ff 'Sage 2.2.1 + ELF10'
    -ff 'Sage 2.2.1+AM1-BCC'
    -ff 'Sage 2.2.1 + AshGC'
    -ff 'Sage 2.3.0rc1 + AshGC'
)
# python plot-aggregated-statistics.py \
#     "${RBFE_FF_KWARGS[@]}" \
#     -i output/jacs/ddG/stats/Sage-2.2.1-+-ELF10.csv \
#     -i output/jacs/ddG/stats/Sage-2.2.1+AM1-BCC.csv \
#     -i output/jacs/ddG/stats/Sage-2.2.1-+-AshGC.csv \
#     -i output/jacs/ddG/stats/Sage-2.3.0rc1-+-AshGC.csv \
#     -im images/compare-refit/jacs-ddG-stats \
#     -ff 'Sage 2.2.1 + ELF10' -ff 'Sage 2.2.1+AM1-BCC' -ff 'Sage 2.2.1 + AshGC' -ff 'Sage 2.3.0rc1 + AshGC' \
#     -u '[kcal/mol]' > logs/plot/aggregated-stats_jacs-ddG-refit.log

# python plot-aggregated-statistics.py \
#     "${RBFE_FF_KWARGS[@]}" \
#     -i output/jacs/dG/stats/Sage-2.2.1-+-ELF10.csv \
#     -i output/jacs/dG/stats/Sage-2.2.1+AM1-BCC.csv \
#     -i output/jacs/dG/stats/Sage-2.2.1-+-AshGC.csv \
#     -i output/jacs/dG/stats/Sage-2.3.0rc1-+-AshGC.csv \
#     -im images/compare-refit/jacs-dG-stats \
#     -ff 'Sage 2.2.1 + ELF10' -ff 'Sage 2.2.1+AM1-BCC' -ff 'Sage 2.2.1 + AshGC' -ff 'Sage 2.3.0rc1 + AshGC' \
#     -u '[kcal/mol]' > logs/plot/aggregated-stats_jacs-dG-refit.log