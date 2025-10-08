#!/usr/bin/env bash

mkdir -p logs/stats

# === qm ===

rm -f logs/stats/qm-rmsds.log

# # === phys-prop ===

# rm -f logs/stats/density-all.log
# rm -f logs/stats/density-validation.log
# rm -f logs/stats/dhmix-all.log
# rm -f logs/stats/dhmix-validation.log

# PHYSPROP_KWARGS=(
#     -cy value
#     -cx reference_value
#     -cye uncertainty
#     -cxe reference_uncertainty
# )

# for FF in 'Sage 2.0.0' \
#           'Sage 2.2.1' \
#           'Sage 2.3.0rc1' \
#           'v1-k100' \
#           'v3-k100' \
#           'smee-v2'
# do
# for FF in 'Sage 2.2.1 + AshGC' ; do

# for FF in 'Sage 2.1.0' \
#           'v0-k20' \
#           'Null 4-mer AAQAA3-2 + OPC3' \
#           'Null AAQAA3-3 + OPC3'
# do
#     datetime=$(date +%Y-%m-%d_%H-%M-%S)
#     echo "[$datetime] Starting statistics computation for phys-prop with force field: $FF"
#     # replace spaces with dashes for filenames
#     FF_FN=$(echo $FF | tr ' ' '-')
#     echo "Force field filename: $FF_FN"
#     python compute-statistics.py                    \
#         "${PHYSPROP_KWARGS[@]}"                      \
#         -ff "$FF"                                   \
#         -i output/phys-prop/density-all/labelled-data.csv   \
#         -o output/phys-prop/density-all/stats/$FF_FN.csv    \
#         >> logs/stats/density-all.log

#     python compute-statistics.py                    \
#         "${PHYSPROP_KWARGS[@]}"                      \
#         -ff "$FF"                                   \
#         -i output/phys-prop/density-validation/labelled-data.csv   \
#         -o output/phys-prop/density-validation/stats/$FF_FN.csv    \
#         >> logs/stats/density-validation.log

#     python compute-statistics.py                    \
#         "${PHYSPROP_KWARGS[@]}"                      \
#         -ff "$FF"                                   \
#         -i output/phys-prop/dhmix-all/labelled-data.csv   \
#         -o output/phys-prop/dhmix-all/stats/$FF_FN.csv    \
#         >> logs/stats/dhmix-all.log

#     python compute-statistics.py                    \
#         "${PHYSPROP_KWARGS[@]}"                      \
#         -ff "$FF"                                   \
#         -i output/phys-prop/dhmix-validation/labelled-data.csv   \
#         -o output/phys-prop/dhmix-validation/stats/$FF_FN.csv    \
#         >> logs/stats/dhmix-validation.log

# done


# # === sfes ===
SFE_KWARGS=(
    -cy 'Value (kcal / mol)'
    -cx 'Reference Value (kcal / mol)'
    -cye 'Uncertainty (kcal / mol)'
    -cxe 'Reference Uncertainty (kcal / mol)'
)

rm -f logs/stats/freesolv.log
rm -f logs/stats/mnsol.log
rm -f logs/stats/tfes.log

# for FF in 'Sage 2.2.1 + AmberTools' \
#           'Sage 2.2.1 + AshGC' \
#           'Sage 2.2.1 + ELF10' \
#           'Sage 2.3.0rc1' \
#           'v1-k100' \
#           'v3-k100' \
#           'v0-k20' ;
# do
for FF in  "fit-iter-1" "vdw-refit" ; do
    datetime=$(date +%Y-%m-%d_%H-%M-%S)
    echo "[$datetime] Starting statistics computation for SFEs with force field: $FF"
    # replace spaces with dashes for filenames
    FF_FN=$(echo $FF | tr ' ' '-')
    echo "Force field filename: $FF_FN"

    python compute-statistics.py                    \
        "${SFE_KWARGS[@]}"                          \
        -ff "$FF"                                   \
        -i output/sfes/mnsol/labelled-data.csv      \
        -o output/sfes/mnsol/stats/$FF_FN.csv       \
        >> logs/stats/mnsol.log

    python compute-statistics.py                    \
        "${SFE_KWARGS[@]}"                          \
        -ff "$FF"                                   \
        -i output/sfes/freesolv/labelled-data.csv   \
        -o output/sfes/freesolv/stats/$FF_FN.csv    \
        >> logs/stats/freesolv.log
    
    python compute-statistics.py                    \
        "${SFE_KWARGS[@]}"                          \
        -ff "$FF"                                   \
        -i output/sfes/tfes/labelled-data.csv       \
        -o output/sfes/tfes/stats/$FF_FN.csv        \
        >> logs/stats/tfes.log

    datetime=$(date +%Y-%m-%d_%H-%M-%S)
    echo "[$datetime] Completed statistics computation for SFEs with force field: $FF".
done



# # === jacs ===
# JACS_KWARGS=(
#     -cy 'Value (kcal / mol)'
#     -cx 'Reference Value (kcal / mol)'
#     -cye 'Uncertainty (kcal / mol)'
#     -cxe 'Reference Uncertainty (kcal / mol)'
# )
# rm -f logs/stats/rbfes-dG.log
# rm -f logs/stats/rbfes-ddG.log
# for FF in 'Sage 2.2.1 + ELF10' \
#           'Sage 2.2.1 + AshGC' \
#           'Sage 2.3.0rc1 + AshGC' \
#           'Sage 2.2.1+AM1-BCC' ;
# do
#     datetime=$(date +%Y-%m-%d_%H-%M-%S)
#     echo "[$datetime] Starting statistics computation for JACS with force field: $FF"
#     # replace spaces with dashes for filenames
#     FF_FN=$(echo $FF | tr ' ' '-')
#     echo "Force field filename: $FF_FN"
#     python compute-statistics.py                    \
#         "${JACS_KWARGS[@]}"                          \
#         -ff "$FF"                                   \
#         -i output/jacs/dG/labelled-properties.csv         \
#         -o output/jacs/dG/stats/$FF_FN.csv          \
#         >> logs/stats/rbfes-dG.log
    
#     python compute-statistics.py                    \
#         "${JACS_KWARGS[@]}"                          \
#         -ff "$FF"                                   \
#         -i output/jacs/ddG/labelled-properties.csv        \
#         -o output/jacs/ddG/stats/$FF_FN.csv         \
#         >> logs/stats/rbfes-ddG.log

# done