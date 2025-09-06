#!/usr/bin/env bash

FORCEFIELD_DIR=../04_benchmark/forcefields

mkdir -p logs/ratio

python get-overrepresented-groups.py \
    -i  ../04_benchmark/phys-prop/output/validation/summary-benchmarks-Density.csv  \
    -nf 'Sage 2.0.0' 'openff-2.0.0'   \
    -nf 'Sage 2.2.1' 'openff-2.2.1'   \
    -nf 'Sage 2.3.0rc1' 'openff-2.3.0rc1'   \
    -nf 'v1-k100'    'fb-fit-v1-single-mean-k100'   \
    -nf 'v3-k100'    'fb-fit-v3-single-mean-k100'   \
    -nf 'smee-v2'   'smee-spice2-systematic-torsion-generation_constrained-linear'  \
    -o  output/phys-prop/density-validation    \
    -s  smiles_1    -s  smiles_2    \
    -e  0.05        -m  abs         \
    -cv 'value'     -crv 'reference_value'         \
    -cm 'forcefield'    -cff FF         \
    -cid 'id'  \
    -lcm    labels/checkmol                 \
    -lff    labels/forcefields/v1-k100      \
    -ff     "${FORCEFIELD_DIR}/fb-fit-v1-single-mean-k100.offxml"      \
    --restrict-to-shared-observables > logs/ratio/density-validation.log

python get-overrepresented-groups.py \
    -i  ../04_benchmark/phys-prop/output/validation/summary-benchmarks-EnthalpyOfMixing.csv  \
    -nf 'Sage 2.0.0' 'openff-2.0.0'   \
    -nf 'Sage 2.2.1' 'openff-2.2.1'   \
    -nf 'Sage 2.3.0rc1' 'openff-2.3.0rc1'   \
    -nf 'v1-k100'    'fb-fit-v1-single-mean-k100'   \
    -nf 'v3-k100'    'fb-fit-v3-single-mean-k100'   \
    -nf 'smee-v2'   'smee-spice2-systematic-torsion-generation_constrained-linear'  \
    -o  output/phys-prop/dhmix-validation    \
    -s  smiles_1    -s  smiles_2    \
    -e  0.3        -m  abs         \
    -cv 'value'     -crv 'reference_value'         \
    -cm 'forcefield'    -cff FF         \
    -cid 'id'  \
    -lcm    labels/checkmol                 \
    -lff    labels/forcefields/v1-k100      \
    -ff     "${FORCEFIELD_DIR}/fb-fit-v1-single-mean-k100.offxml"      \
    --restrict-to-shared-observables > logs/ratio/dhmix-validation.log
    

python get-overrepresented-groups.py \
    -i  ../04_benchmark/phys-prop/output/training/summary-benchmarks-Density.csv    \
    -i  ../04_benchmark/phys-prop/output/validation/summary-benchmarks-Density.csv  \
    -nf 'Sage 2.0.0' 'openff-2.0.0'   \
    -nf 'Sage 2.2.1' 'openff-2.2.1'   \
    -nf 'Sage 2.3.0rc1' 'openff-2.3.0rc1'   \
    -nf 'v1-k100'    'fb-fit-v1-single-mean-k100'   \
    -nf 'v3-k100'    'fb-fit-v3-single-mean-k100'   \
    -nf 'smee-v2'   'smee-spice2-systematic-torsion-generation_constrained-linear'  \
    -o  output/phys-prop/density    \
    -s  smiles_1    -s  smiles_2    \
    -e  0.05        -m  abs         \
    -cv 'value'     -crv 'reference_value'         \
    -cm 'forcefield'    -cff FF         \
    -cid 'id'  \
    -lcm    labels/checkmol                 \
    -lff    labels/forcefields/v1-k100      \
    -ff     "${FORCEFIELD_DIR}/fb-fit-v1-single-mean-k100.offxml"      \
    --restrict-to-shared-observables > logs/ratio/density.log

python get-overrepresented-groups.py \
    -i  ../04_benchmark/phys-prop/output/training/summary-benchmarks-EnthalpyOfMixing.csv    \
    -i  ../04_benchmark/phys-prop/output/validation/summary-benchmarks-EnthalpyOfMixing.csv  \
    -nf 'Sage 2.0.0' 'openff-2.0.0'   \
    -nf 'Sage 2.2.1' 'openff-2.2.1'   \
    -nf 'Sage 2.3.0rc1' 'openff-2.3.0rc1'   \
    -nf 'v1-k100'    'fb-fit-v1-single-mean-k100'   \
    -nf 'v3-k100'    'fb-fit-v3-single-mean-k100'   \
    -nf 'smee-v2'   'smee-spice2-systematic-torsion-generation_constrained-linear'  \
    -o  output/phys-prop/dhmix    \
    -s  smiles_1    -s  smiles_2    \
    -e  0.3        -m  abs         \
    -cv 'value'     -crv 'reference_value'         \
    -cm 'forcefield'    -cff FF         \
    -cid 'id'  \
    -lcm    labels/checkmol                 \
    -lff    labels/forcefields/v1-k100      \
    -ff     "${FORCEFIELD_DIR}/fb-fit-v1-single-mean-k100.offxml"      \
    --restrict-to-shared-observables > logs/ratio/dhmix.log
    
python get-overrepresented-groups.py \
    -i  ../04_benchmark/rbfes/output/ddG.csv    \
    -nf 'Sage 2.2.1 + ELF10' 'Sage 2.2.1 + ELF10'   \
    -nf 'Sage 2.2.1+AM1-BCC' 'Sage 2.2.1+AM1-BCC'   \
    -nf 'Sage 2.2.1 + AshGC' 'Sage 2.2.1 + AshGC'   \
    -nf 'Sage 2.3.0rc0 + ELF10' 'Sage 2.3.0rc0 + ELF10' \
    -nf 'Sage 2.3.0rc0 + AshGC' 'Sage 2.3.0rc0 + AshGC' \
    -nf 'Sage 2.3.0rc1 + AshGC' 'Sage 2.3.0rc1 + AshGC' \
    -o  output/jacs/ddG    \
    -s  'SMILES 1'    -s  'SMILES 2'    \
    -e  1.2       -m  abs         \
    -cv 'Value (kcal / mol)'     -crv 'Reference Value (kcal / mol)'         \
    -cm 'Force Field'    -cff FF         \
    -cid 'Transformation'  \
    -lcm    labels/checkmol                 \
    -lff    labels/forcefields/v1-k100      \
    -ff     "${FORCEFIELD_DIR}/fb-fit-v1-single-mean-k100.offxml"      \
    --restrict-to-shared-observables > logs/ratio/rbfes-ddG.log
    
python get-overrepresented-groups.py \
    -i  ../04_benchmark/rbfes/output/dG.csv    \
    -nf 'Sage 2.2.1 + ELF10' 'Sage 2.2.1 + ELF10'   \
    -nf 'Sage 2.2.1+AM1-BCC' 'Sage 2.2.1+AM1-BCC'   \
    -nf 'Sage 2.2.1 + AshGC' 'Sage 2.2.1 + AshGC'   \
    -nf 'Sage 2.3.0rc0 + ELF10' 'Sage 2.3.0rc0 + ELF10' \
    -nf 'Sage 2.3.0rc0 + AshGC' 'Sage 2.3.0rc0 + AshGC' \
    -nf 'Sage 2.3.0rc1 + AshGC' 'Sage 2.3.0rc1 + AshGC' \
    -o  output/jacs/dG    \
    -s  'SMILES'    \
    -e  1.2       -m  abs         \
    -cv 'Value (kcal / mol)'     -crv 'Reference Value (kcal / mol)'         \
    -cm 'Force Field'    -cff FF         \
    -cid 'Mapped SMILES'  \
    -lcm    labels/checkmol                 \
    -lff    labels/forcefields/v1-k100      \
    -ff     "${FORCEFIELD_DIR}/fb-fit-v1-single-mean-k100.offxml"      \
    --restrict-to-shared-observables > logs/ratio/rbfes-dG.log

