#!/usr/bin/env bash

mkdir -p logs/plot

REFERENCE="Sage 2.2.1"

# for TARGET in v3-k100 "Sage 2.3.0rc1" "Sage 2.2.1 + AshGC" ; do

#     python plot-change-in-target.py \
#         -i output/phys-prop/density-all/labelled-data.csv \
#         -o "output/phys-prop/density-all/${TARGET}-vs-${REFERENCE}.csv" \
#         -im "images/compare-refit/density-changes/${TARGET}-vs-${REFERENCE}.png" \
#         -ref "$REFERENCE"       \
#         -tgt "$TARGET"          \
#         -cff "FF"  -cv "value" \
#         -cs "smiles_1" -cs "smiles_2" \
#         -pn "Density"  -pu "g/mL" > "logs/plot/density-${TARGET}-vs-${REFERENCE}.log"
# done


for TARGET in v3-k100 "Sage 2.3.0rc1" "Sage 2.2.1 + AshGC" ; do

    python plot-change-in-target.py \
        -i output/phys-prop/dhmix-all/labelled-data.csv \
        -o "output/phys-prop/dhmix-all/${TARGET}-vs-${REFERENCE}.csv" \
        -im "images/compare-refit/dhmix-changes/${TARGET}-vs-${REFERENCE}.png" \
        -ref "$REFERENCE"       \
        -tgt "$TARGET"          \
        -cff "FF"  -cv "value" \
        -cs "smiles_1" -cs "smiles_2" \
        -pn "∆Hmix"  -pu "kJ/mol" > "logs/plot/dhmix-${TARGET}-vs-${REFERENCE}.log"
done

