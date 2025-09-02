#!/usr/bin/env bash

mkdir -p logs

# RBFEs
python label-with-checkmol.py               \
    -i ../04_benchmark/rbfes/output/ddG.csv \
    -s "SMILES 1"                           \
    -s "SMILES 2"                           \
    -o labels/checkmol/rbfes-jacs.parquet > logs/label-checkmol-rbfes.log

python label-with-forcefield.py                 \
    -i   ../04_benchmark/rbfes/output/ddG.csv   \
    -s   "SMILES 1"                             \
    -s   "SMILES 2"                             \
    -ff  ../04_benchmark/forcefields/fb-fit-v1-single-mean-k100.offxml \
    -ffn "v1-k100"                              \
    -o labels/forcefields/v1-k100/rbfes-jacs.parquet > logs/label-ff-v1-k100-rbfes.log

python label-with-forcefield.py                 \
    -i   ../04_benchmark/rbfes/output/ddG.csv   \
    -s   "SMILES 1"                             \
    -s   "SMILES 2"                             \
    -ff  ../04_benchmark/forcefields/fb-fit-v3-single-mean-k100.offxml \
    -ffn "v3-k100"                              \
    -o labels/forcefields/v3-k100/rbfes-jacs.parquet > logs/label-ff-v1-k100-rbfes.log


# SFEs
python label-with-checkmol.py                   \
    -i ../04_benchmark/sfes/output/freesolv.csv \
    -s "Solute"                                 \
    -s "Solvent"                                \
    -o labels/checkmol/sfes-freesolv.parquet > logs/label-checkmol-sfes-freesolv.log

python label-with-forcefield.py                 \
    -i   ../04_benchmark/sfes/output/freesolv.csv   \
    -s   "Solute"                             \
    -s   "Solvent"                             \
    -ff  ../04_benchmark/forcefields/fb-fit-v1-single-mean-k100.offxml \
    -ffn "v1-k100"                              \
    -o labels/forcefields/v1-k100/sfes-freesolv.parquet > logs/label-ff-v1-k100-sfes-freesolv.log

python label-with-forcefield.py                 \
    -i   ../04_benchmark/sfes/output/freesolv.csv   \
    -s   "Solute"                             \
    -s   "Solvent"                             \
    -ff  ../04_benchmark/forcefields/fb-fit-v3-single-mean-k100.offxml \
    -ffn "v3-k100"                              \
    -o labels/forcefields/v3-k100/sfes-freesolv.parquet > logs/label-ff-v3-k100-sfes-freesolv.log


python label-with-checkmol.py                   \
    -i ../04_benchmark/sfes/output/mnsol.csv \
    -s "Solute"                                 \
    -s "Solvent"                                \
    -o labels/checkmol/sfes-mnsol.parquet > logs/label-checkmol-sfes-mnsol.log

python label-with-forcefield.py                 \
    -i   ../04_benchmark/sfes/output/mnsol.csv   \
    -s   "Solute"                             \
    -s   "Solvent"                             \
    -ff  ../04_benchmark/forcefields/fb-fit-v1-single-mean-k100.offxml \
    -ffn "v1-k100"                              \
    -o labels/forcefields/v1-k100/sfes-mnsol.parquet > logs/label-ff-v1-k100-sfes-mnsol.log

python label-with-forcefield.py                 \
    -i   ../04_benchmark/sfes/output/mnsol.csv   \
    -s   "Solute"                             \
    -s   "Solvent"                             \
    -ff  ../04_benchmark/forcefields/fb-fit-v3-single-mean-k100.offxml \
    -ffn "v3-k100"                              \
    -o labels/forcefields/v3-k100/sfes-mnsol.parquet > logs/label-ff-v3-k100-sfes-mnsol.log

