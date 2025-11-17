#!/usr/bin/env bash

mkdir logs
# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solvent/v3-vs-221+elf10"     \
#     -r  "Sage 2.2.1 + ELF10"                        \
#     -t  "v3-k100"                                   \
#     -cp Solvent                                     \
#     -co Solute                                      \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solvent-v3-vs-221+elf10.log

# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solute/v3-vs-221+elf10"      \
#     -r  "Sage 2.2.1 + ELF10"                        \
#     -t  "v3-k100"                                   \
#     -cp Solute                                      \
#     -co Solvent                                     \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solute-v3-vs-221+elf10.log

# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solvent/v3-vs-221+ashgc"   \
#     -r  "Sage 2.2.1 + AshGC"                        \
#     -t  "v3-k100"                                   \
#     -cp Solvent                                     \
#     -co Solute                                      \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solvent-v3-vs-221+ashgc.log

# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solute/v3-vs-221+ashgc"    \
#     -r  "Sage 2.2.1 + AshGC"                        \
#     -t  "v3-k100"                                   \
#     -cp Solute                                      \
#     -co Solvent                                     \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solute-v3-vs-221+ashgc.log


# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solvent/221+ashgc-vs-221+elf10"   \
#     -r  "Sage 2.2.1 + ELF10"                        \
#     -t  "Sage 2.2.1 + AshGC"                        \
#     -cp Solvent                                     \
#     -co Solute                                      \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solvent-221-ashgc-vs-elf10.log

# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solute/221+ashgc-vs-221+elf10"    \
#     -r  "Sage 2.2.1 + ELF10"                        \
#     -t  "Sage 2.2.1 + AshGC"                        \
#     -cp Solute                                      \
#     -co Solvent                                     \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solute-221-ashgc-vs-elf10.log


# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solute/rc1-vs-221+elf10"    \
#     -r  "Sage 2.2.1 + ELF10"                        \
#     -t  "Sage 2.3.0rc1"                        \
#     -cp Solute                                      \
#     -co Solvent                                     \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solute-rc1-vs-221-elf10.log

# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solute/v3-vs-rc1"    \
#     -t  "v3-k100"                        \
#     -r  "Sage 2.3.0rc1"                        \
#     -cp Solute                                      \
#     -co Solvent                                     \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solute-v3-vs-rc1.log


# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solute/v3-vs-v0"    \
#     -t  "v3-k100"                        \
#     -r  "v0-k20"                        \
#     -cp Solute                                      \
#     -co Solvent                                     \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solute-v3-vs-v0.log

# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solute/v0-vs-221+elf10"    \
#     -r  "Sage 2.2.1 + ELF10"                        \
#     -t  "v0-k20"                        \
#     -cp Solute                                      \
#     -co Solvent                                     \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solute-v0-vs-221-elf10.log


# python plot-difference-graph.py                     \
#     -i  ../../../04_benchmark/sfes/output/mnsol.csv \
#     -o  "images/mnsol-by-solute/iter-1-vs-221+elf10"    \
#     -r  "Sage 2.2.1 + ELF10"                        \
#     -t  "fit-iter-1"                        \
#     -cp Solute                                      \
#     -co Solvent                                     \
#     -w  6                                           \
#     -m  50                                          \
#     -x  5 > logs/mnsol-by-solute-iter-1-vs-221-elf10.log

python plot-difference-graph.py                     \
    -i  ../../../04_benchmark/sfes/output/mnsol.csv \
    -o  "images/mnsol-by-solute/iter-1-vs-ashgc"    \
    -r  "Sage 2.2.1 + AshGC"                        \
    -t  "fit-iter-1"                        \
    -cp Solute                                      \
    -co Solvent                                     \
    -w  6                                           \
    -m  50                                          \
    -x  5 > logs/mnsol-by-solute-iter-1-vs-ashgc.log
