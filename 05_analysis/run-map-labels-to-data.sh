#!/usr/bin/env bash

mkdir -p logs/map

# combine and label data computed in `04_benchmark`

# === physical properties ===
# PHYSPROP_DIR="../04_benchmark/phys-prop/output"
# PHYSPROP_FFS=(
#     -nf 'Sage 2.0.0' 'openff-2.0.0'
#     -nf 'Sage 2.1.0' 'openff-2.1.0'
#     -nf 'Sage 2.2.1' 'openff-2.2.1'
#     -nf 'Sage 2.3.0rc1' 'openff-2.3.0rc1'
#     -nf 'v0-k20'     'fb-fit-v0-single-mean-k20'
#     -nf 'v1-k100'    'fb-fit-v1-single-mean-k100'
#     -nf 'v3-k100'    'fb-fit-v3-single-mean-k100'
#     -nf 'smee-v2'   'smee-spice2-systematic-torsion-generation_constrained-linear'
#     -nf 'Sage 2.2.1 + AshGC' 'openff-2.2.1-ashgc'
#     -nf 'Null 4-mer AAQAA3-2 + OPC3' 'Protein-Null-0.0.3-4-mer-AAQAA3-1E3-2-OPC3_opc3'
#     -nf 'Null AAQAA3-3 + OPC3' 'Protein-Null-0.0.3-Pair-NMR-1E4-Umbrella-1E3-4-AAQAA3-1E3-3-OPC3_opc3'
# )
# PHYSPROP_KWARGS=(
#     -s  smiles_1    -s  smiles_2
#     -cm 'forcefield'    -cff FF
# )
    
# python map-labels-to-data.py        \
#     "${PHYSPROP_FFS[@]}" "${PHYSPROP_KWARGS[@]}"  \
#     -i  ${PHYSPROP_DIR}/training/summary-benchmarks-Density.csv     \
#     -i  ${PHYSPROP_DIR}/validation/summary-benchmarks-Density.csv   \
#     -o  output/phys-prop/density-all/labelled-data.csv              \
#     > logs/map/map-phys_prop-density_all.log

# python map-labels-to-data.py        \
#     "${PHYSPROP_FFS[@]}" "${PHYSPROP_KWARGS[@]}"  \
#     -i  ${PHYSPROP_DIR}/training/summary-benchmarks-EnthalpyOfMixing.csv    \
#     -i  ${PHYSPROP_DIR}/validation/summary-benchmarks-EnthalpyOfMixing.csv  \
#     -o  output/phys-prop/dhmix-all/labelled-data.csv                        \
#     > logs/map/map-phys_prop-dhmix_all.log

# python map-labels-to-data.py        \
#     "${PHYSPROP_FFS[@]}" "${PHYSPROP_KWARGS[@]}"  \
#     -i  ${PHYSPROP_DIR}/validation/summary-benchmarks-Density.csv   \
#     -o  output/phys-prop/density-validation/labelled-data.csv       \
#     > logs/map/map-phys_prop-density_validation.log

# python map-labels-to-data.py        \
#     "${PHYSPROP_FFS[@]}" "${PHYSPROP_KWARGS[@]}"  \
#     -i  ${PHYSPROP_DIR}/validation/summary-benchmarks-EnthalpyOfMixing.csv  \
#     -o  output/phys-prop/dhmix-validation/labelled-data.csv                 \
#     > logs/map/map-phys_prop-dhmix_validation.log


# # === sfes ===

# SFE_FFS=(
#     -nf 'Sage 2.2.1 + ELF10' 'Sage 2.2.1 + ELF10'
#     -nf 'Sage 2.2.1 + AmberTools' 'Sage 2.2.1 + AmberTools'
#     -nf 'Sage 2.2.1 + AshGC' 'Sage 2.2.1 + AshGC'
#     -nf 'Sage 2.3.0rc1' 'Sage 2.3.0rc1'
#     -nf 'fit-iter-1' 'fit-iter-1'
#     -nf 'vdw-refit' 'vdw-refit'
#     -nf 'v0-k20' 'v0-k20'
#     -nf 'v1-k100' 'v1-k100'
#     -nf 'v3-k100' 'v3-k100'
# )
# SFE_KWARGS=(
#     -cm 'Method'    -cff FF
# )


# python map-labels-to-data.py                        \
#     "${SFE_FFS[@]}" "${SFE_KWARGS[@]}"              \
#     -s  'Solute'  -s 'Solvent' \
#     -i  ../04_benchmark/sfes/output/freesolv.csv    \
#     -o  output/sfes/freesolv/labelled-data.csv      \
#     > logs/map/map-sfes-freesolv.log

# python map-labels-to-data.py                        \
#     "${SFE_FFS[@]}" "${SFE_KWARGS[@]}"              \
#     -s  'Solute'  -s 'Solvent' \
#     -i  ../04_benchmark/sfes/output/mnsol.csv       \
#     -o  output/sfes/mnsol/labelled-data.csv         \
#     > logs/map/map-sfes-mnsol.log

# python map-labels-to-data.py                        \
#     "${SFE_FFS[@]}" "${SFE_KWARGS[@]}"              \
#     -s  'Solute'  -s 'Solvent 2' \
#     -i  ../04_benchmark/sfes/output/tfes.csv        \
#     -o  output/sfes/tfes/labelled-data.csv          \
#     > logs/map/map-sfes-tfes.log


# sfe solutes

# SFE_FFS=(
#     -nf 'Sage 2.2.1 + ELF10' 'Sage 2.2.1 + ELF10'
#     -nf 'Sage 2.2.1 + AmberTools' 'Sage 2.2.1 + AmberTools'
#     -nf 'Sage 2.2.1 + AshGC' 'Sage 2.2.1 + AshGC'
#     -nf 'Sage 2.3.0rc1' 'Sage 2.3.0rc1'
#     -nf 'fit-iter-1' 'fit-iter-1'
#     -nf 'vdw-refit' 'vdw-refit'
#     -nf 'v0-k20' 'v0-k20'
#     -nf 'v1-k100' 'v1-k100'
#     -nf 'v3-k100' 'v3-k100'
# )
# SFE_KWARGS=(
#     -cm 'Method'    -cff FF
# )


# python map-labels-to-data.py                        \
#     "${SFE_FFS[@]}" "${SFE_KWARGS[@]}"              \
#     -s  'Solute'  \
#     -i  ../04_benchmark/sfes/output/freesolv.csv    \
#     -o  output/sfes-solutes/freesolv/labelled-data.csv      \
#     > logs/map/map-sfes-solutes-freesolv.log

# python map-labels-to-data.py                        \
#     "${SFE_FFS[@]}" "${SFE_KWARGS[@]}"              \
#     -s  'Solute'  \
#     -i  ../04_benchmark/sfes/output/mnsol.csv       \
#     -o  output/sfes-solutes/mnsol/labelled-data.csv         \
#     > logs/map/map-sfes-solutes-mnsol.log

# python map-labels-to-data.py                        \
#     "${SFE_FFS[@]}" "${SFE_KWARGS[@]}"              \
#     -s  'Solute'  \
#     -i  ../04_benchmark/sfes/output/tfes.csv        \
#     -o  output/sfes-solutes/tfes/labelled-data.csv          \
#     > logs/map/map-sfes-solutes-tfes.log
# === RBFEs ===

RBFE_FFS=(
    -nf 'Sage 2.2.1 + ELF10' 'Sage 2.2.1 + ELF10'
    -nf 'Sage 2.2.1+AM1-BCC' 'Sage 2.2.1 + AM1-BCC'
    -nf 'Sage 2.2.1 + AshGC' 'Sage 2.2.1 + AshGC'
    -nf 'Sage 2.3.0rc0 + ELF10' 'Sage 2.3.0rc0 + ELF10'
    -nf 'Sage 2.3.0rc0 + AshGC' 'Sage 2.3.0rc0 + AshGC'
    -nf 'Sage 2.3.0rc1 + ELF10' 'Sage 2.3.0rc1 + ELF10'
    -nf 'Sage 2.3.0rc1 + AshGC' 'Sage 2.3.0rc1 + AshGC'
    -nf 'Sage 2.3.0rc2 + AshGC' 'Sage 2.3.0rc2 + AshGC'
    -nf 'Sage 2.3.0rc2 + ELF10' 'Sage 2.3.0rc2 + ELF10'
)
RBFE_KWARGS=(
    -cm 'Force Field'    -cff FF
)

python map-labels-to-data.py                        \
    "${RBFE_FFS[@]}" "${RBFE_KWARGS[@]}"            \
    -s 'SMILES'   \
    -i  ../04_benchmark/rbfes/output/dG.csv         \
    -o  output/jacs/dG/labelled-data.csv           \
    > logs/map/map-rbfes-dG.log

python map-labels-to-data.py                        \
    "${RBFE_FFS[@]}" "${RBFE_KWARGS[@]}"            \
    -s 'SMILES 1' -s 'SMILES 2'   \
    -i  ../04_benchmark/rbfes/output/ddG.csv         \
    -o  output/jacs/ddG/labelled-data.csv           \
    > logs/map/map-rbfes-ddG.log

# === QM ===

# QM_FFS=(
#     -nf 'Sage 2.0.0' 'openff_unconstrained-2.0.0'
#     -nf 'Sage 2.2.1' 'openff_unconstrained-2.2.1'
#     -nf 'Sage 2.3.0rc1' 'openff_unconstrained-2.3.0rc1'
#     -nf 'v1-k100'    'fb-fit-v1-single-mean-k100_unconstrained'
#     -nf 'v3-k100'    'fb-fit-v3-single-mean-k100_unconstrained'
# )
# QM_KWARGS=(
#     -cm 'method'    -cff FF  -s 'mapped_smiles'
# )

# python map-labels-to-data.py                        \
#     "${QM_FFS[@]}" "${QM_KWARGS[@]}"                \
#     -i  ../04_benchmark/qm/all-to-all-rmsd-tfd      \
#     -o  output/qm/all-to-all-rmsd-tfd/labelled-data.csv \
#     > logs/map/map-qm-all-to-all-rmsd-tfd.log

# python map-labels-to-data.py                        \
#     "${QM_FFS[@]}" "${QM_KWARGS[@]}"                \
#     -i  ../04_benchmark/qm/ddes \
#     -o  output/qm/ddes/labelled-data.csv            \
#     > logs/map/map-qm-ddes.log



# python map-labels-to-data.py                        \
#     "${SFE_FFS[@]}" "${SFE_KWARGS[@]}"              \
#     -s  'smiles_1'  -s 'smiles_2' \
#     -i  ../02_fit-vdw/investigate-refit/output/sage-gradient-contributions.csv    \
#     -o  ../02_fit-vdw/investigate-refit/output/sage-gradient-contributions-labelled.csv