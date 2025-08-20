#!/bin/bash
mkdir -p logs output counts

DATA_DIR="../../01_download-data/qm/data/tables"

# code below curates data for particular versions of valence force fields
# it selects torsiondrives and optimizations from the QM data
# two sets of optimizations are selected -- one with only a single lowest-energy conformer
# per molecule, and one with multiple conformers per molecule.

# In experimentations with 2.3.0rc1, we found that the single-conformer datasets typically
# performed slightly better on QM benchmarks.

for VER in v0 v1 v2 v3; do

    echo "==== Selecting for ${VER} ===="

    python select-torsiondrives-smallest.py  \
        -t  "${DATA_DIR}/torsiondrive" \
        -xc bad-data/failed-charge-cmiles.dat    \
        -xs bad-data/smarts-to-exclude.dat \
        -xq bad-data/failed-seesaw-sx4-tds.dat \
        -i  "label/parameters/torsions/ash-sage-valence-${VER}" \
        -ff "../01_generate-forcefield/output/initial-force-field-${VER}.offxml" \
        -o  "output/torsiondrives-${VER}.json" \
        -n 12 -np 8 \
        -oc "counts/torsion-counts-${VER}.json" > logs/select-torsiondrives-${VER}.log 2>&1

    python select-optimizations-smallest.py  \
        -t "${DATA_DIR}/optimization" \
        -xs bad-data/smarts-to-exclude.dat \
        -xq bad-data/failed-seesaw-sx4.dat \
        -xq  bad-data/bad-qcarchive_ids.dat           \
        -xq  bad-data/existing-filtered-ids.dat       \
        -xc bad-data/failed-charge-cmiles.dat    \
        -xdn "OpenFF Industry Benchmark Season 1 v1.2" \
        -i  "label/parameters/valence/ash-sage-valence-${VER}" \
        -ff "../01_generate-forcefield/output/initial-force-field-${VER}.offxml" \
        -o  "output/optimizations-single-${VER}.json" \
        -nc 1  -n 100 -np 8 \
        -oc "counts/valence-counts-single-${VER}.json" > logs/select-optimizations-single-${VER}.log 2>&1


    python select-optimizations-smallest.py  \
        -t "${DATA_DIR}/optimization" \
        -xs bad-data/smarts-to-exclude.dat \
        -xq bad-data/failed-seesaw-sx4.dat \
        -xq  bad-data/bad-qcarchive_ids.dat           \
        -xq  bad-data/existing-filtered-ids.dat       \
        -xc bad-data/failed-charge-cmiles.dat    \
        -xdn "OpenFF Industry Benchmark Season 1 v1.2" \
        -i  "label/parameters/valence/ash-sage-valence-${VER}" \
        -ff "../01_generate-forcefield/output/initial-force-field-${VER}.offxml" \
        -o  "output/optimizations-multi-${VER}.json" \
        -nc 10  -n 30 -np 8 \
        -oc "counts/valence-counts-multi-${VER}.json" > logs/select-optimizations-multi-${VER}.log 2>&1

done

# the below code selects validation sets for testing
# this is just shown as an example of how to use the `-xd` flag

python select-torsiondrives-smallest.py             \
    -t  "${DATA_DIR}/torsiondrive"                  \
    -xc bad-data/failed-charge-cmiles.dat           \
    -xs bad-data/smarts-to-exclude.dat              \
    -xq bad-data/failed-seesaw-sx4-tds.dat          \
    -xd output/torsiondrives-v0.json                \
    -xd output/torsiondrives-v1.json                \
    -xd output/torsiondrives-v2.json                \
    -xd output/torsiondrives-v3.json                \
    -i  "label/parameters/torsions/ash-sage-valence-v3"                     \
    -ff "../01_generate-forcefield/output/initial-force-field-v3.offxml"    \
    -o  "output/torsiondrives-validation.json" \
    -n 12 -np 8 \
    -oc "counts/torsion-counts-validation.json" > logs/select-torsiondrives-validation.log 2>&1

python select-optimizations-smallest.py             \
    -t "${DATA_DIR}/optimization"                   \
    -xs bad-data/smarts-to-exclude.dat              \
    -xq bad-data/failed-seesaw-sx4.dat              \
    -xq  bad-data/bad-qcarchive_ids.dat           \
    -xq  bad-data/existing-filtered-ids.dat       \
    -xc bad-data/failed-charge-cmiles.dat           \
    -xd output/optimizations-single-v0.json         \
    -xd output/optimizations-single-v1.json         \
    -xd output/optimizations-single-v2.json         \
    -xd output/optimizations-single-v3.json         \
    -xd output/optimizations-multi-v0.json          \
    -xd output/optimizations-multi-v1.json          \
    -xd output/optimizations-multi-v2.json          \
    -xd output/optimizations-multi-v3.json          \
    -xdn "OpenFF Industry Benchmark Season 1 v1.2" \
    -i  "label/parameters/valence/ash-sage-valence-v3" \
    -ff "../01_generate-forcefield/output/initial-force-field-v3.offxml" \
    -o  "output/optimizations-validation.json" \
    -nc 20  -n 50 -np 8 \
    -oc "counts/valence-counts-validation.json" > logs/select-optimizations-validation.log 2>&1
