#!/usr/bin/env bash

python compare-training-overlap.py                  \
    -t  ash_training_smiles.smi                     \
    -r  /Users/lily/pydev/IndustryBenchmarks2024/industry_benchmarks/input_structures/prepared_structures    \
    -o  output                                      \
    -im images                                      \
    --plot-cores