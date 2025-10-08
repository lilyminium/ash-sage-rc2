#!/usr/bin/env bash

mkdir -p logs/plot

FFDIR="../04_benchmark/forcefields"

python plot-lj-curves.py                                        \
    -rf "openff-2.2.1.offxml"                                   \
    -rn "Sage 2.2.1"                                            \
    -tf "${FFDIR}/fb-fit-v3-single-mean-k100.offxml"            \
    -tn "v3-k100"                                               \
    -id images/compare-refit/lj-curves/sage-2.2.1_vs_v3-k100    \
    > logs/plot/plot-lj-curves_sage-2.2.1_vs_v3-k100.log
