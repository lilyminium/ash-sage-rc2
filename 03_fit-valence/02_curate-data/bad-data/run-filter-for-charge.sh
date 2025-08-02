#!/usr/bin/env bash
#SBATCH -J filter-charge
#SBATCH -p standard
#SBATCH -t 96:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A.out

source ~/.bashrc

conda activate evaluator-050

DIRECTORY="../../../01_download-data/qm/data/tables"

python filter-out-smiles-by-charge.py   \
    -np 64 \
    -i $DIRECTORY/optimization \
    -i $DIRECTORY/torsiondrive \
    -o failed-charge-cmiles.dat
