#!/bin/bash
#SBATCH -J check-properties
#SBATCH -p standard
#SBATCH -t 01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mem=4gb
#SBATCH --constraint=fastscratch
#SBATCH --output slurm-%x.%A.out



. ~/.bashrc

# Use the right conda environment
conda activate evaluator-050

TIER="validation"
TIER="training"

DATASET="../../01_download-data/physprop/final/output/${TIER}-set.json"


python check-properties-to-run.py      \
    -i $DATASET                 \
    -d $TIER                    \
    -o output/$TIER
