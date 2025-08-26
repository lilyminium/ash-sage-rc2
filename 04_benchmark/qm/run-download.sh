#!/usr/bin/env bash
#SBATCH -J download
#SBATCH -p standard
#SBATCH -t 16:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate nagl-valence

DATA_DIR="../../03_fit-valence/02_curate-data/output"


python get-optimization-data.py -i $DATA_DIR/optimizations-validation.json -o data/optimization

python get-torsiondrive-data.py -i $DATA_DIR/torsiondrives-validation.json -o data/torsiondrive

python download-industry-set.py -o data/optimization

