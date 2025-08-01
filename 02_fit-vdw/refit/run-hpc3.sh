#!/usr/bin/env bash
#SBATCH -J slurm-nagl-v2
#SBATCH -p standard
#SBATCH -t 7-00:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --output slurm-%x.%A.out

PORT=8009

source ~/.bashrc

# ENVFILE="evaluator-test-env-openff"
ENVFILE="evaluator-050"

conda activate $ENVFILE

# write force field
# python set-up-forcefield.py -n 3


# run fit
python execute-fit-slurm-distributed.py                 \
    --port                  $PORT                       \
    --n-min-workers         1                           \
    --n-max-workers         23                          \
    --memory-per-worker     4                           \
    --walltime              "08:00:00"                  \
    --queue                 "free-gpu"                       \
    --conda-env             $ENVFILE                    \
    --extra-script-option   "--gres=gpu:1"              # note: this is hpc3 specific
