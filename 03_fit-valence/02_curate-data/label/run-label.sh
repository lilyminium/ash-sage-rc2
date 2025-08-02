#!/bin/bash

echo "=== Labelling === "

FF_DIR="../../01_generate-forcefield/output"
DATA_DIR="../../../01_download-data/qm/data/tables"

mkdir -p logs

for VER in v0 v1 v2 v3; do

    echo "==== Labeling for ${VER} ===="

    python label-valence-with-forcefield.py                                     \
        -i      "${DATA_DIR}/optimization"                                      \
        -ff     "${FF_DIR}/initial-force-field-${VER}.offxml"                   \
        -ffn    "ash-sage-valence-${VER}" > logs/label-valence-${VER}.log 2>&1

    python label-torsions-with-forcefield.py                                    \
        -i      "${DATA_DIR}/torsiondrive"                                      \
        -ff     "${FF_DIR}/initial-force-field-${VER}.offxml"                   \
        -ffn    "ash-sage-valence-${VER}" > logs/label-torsion-${VER}.log 2>&1

done

