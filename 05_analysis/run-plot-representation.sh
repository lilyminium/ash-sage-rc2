#!/usr/bin/env bash

mkdir -p logs/plot

for group in "Checkmol" "Bonds" "Angles" "ProperTorsions" "ImproperTorsions" "vdW"
do
    # python plot-overrepresentation.py \
    #     -i output/phys-prop/density/ratio_representation_0.05.csv \
    #     -o images/overrepresentation/phys-prop/densities-$group.png \
    #     -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'v1-k100' -ff 'v3-k100' \
    #     -d "Densities > 0.05 g/mL" \
    #     -g "$group"

    # python plot-overrepresentation.py \
    #     -i output/phys-prop/dhmix/ratio_representation_0.3.csv \
    #     -o images/overrepresentation/phys-prop/dhmix-$group.png \
    #     -ff 'Sage 2.0.0' -ff 'Sage 2.2.1' -ff 'v1-k100' -ff 'v3-k100' \
    #     -d "∆Hmix > 0.3 kJ/mol" \
    #     -g "$group"

    python plot-overrepresentation.py \
        -i output/jacs/ddG/ratio_representation_1.2.csv \
        -o images/overrepresentation/jacs/ddG-$group.png \
        -ff 'Sage 2.2.1 + ELF10' -ff 'Sage 2.2.1+AM1-BCC' -ff 'Sage 2.2.1 + AshGC' -ff 'Sage 2.3.0rc1 + AshGC' \
        -d "∆∆G > 1.2 kcal/mol" \
        -g "$group"

    python plot-overrepresentation.py \
        -i output/jacs/dG/ratio_representation_10.0.csv \
        -o images/overrepresentation/jacs/dG-$group.png \
        -ff 'Sage 2.2.1 + ELF10' -ff 'Sage 2.2.1+AM1-BCC' -ff 'Sage 2.2.1 + AshGC' -ff 'Sage 2.3.0rc1 + AshGC' \
        -d "∆G > 10 kcal/mol" \
        -g "$group"
done


