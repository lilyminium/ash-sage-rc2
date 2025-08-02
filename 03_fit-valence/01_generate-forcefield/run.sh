#!/bin/bash

mkdir output


VDW_FF="../../02_fit-vdw/refit/result/optimize/force-field.offxml"

# VDW_FF="openff-2.2.1.offxml"

python generate-forcefield-v0.py -i $VDW_FF

python generate-forcefield-v1.py -i $VDW_FF

python generate-forcefield-v2.py -i $VDW_FF

python generate-forcefield-v3.py -i $VDW_FF
