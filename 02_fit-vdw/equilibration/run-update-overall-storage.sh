#!/usr/bin/env bash
#SBATCH -J update-storage
#SBATCH -p standard
#SBATCH --account dmobley_lab
#SBATCH -t 01:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --output slurm-%x.%A.out

. ~/.bashrc

# Use the right conda environment
conda activate evaluator-050


python update-overall-storage.py -lfs stored_data -lsed ../equilibration_stored_data
