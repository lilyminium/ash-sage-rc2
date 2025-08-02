#!/bin/bash
#SBATCH -J export
#SBATCH -p standard
#SBATCH -t 08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mem=16gb
#SBATCH --constraint=fastscratch
#SBATCH --output slurm-%x.%A.out

date
hostname

source ~/.bashrc
conda activate nagl-valence    

OUTDIR="../../04_benchmark/forcefields"

mkdir -p ${OUTDIR}


INPUTFILE="${FFNAME}/result/optimize/force-field.offxml"
OUTPUTFILE="${OUTDIR}/${FFNAME}.offxml"

python remove_cosmetic_attributes.py -i ${INPUTFILE} -o ${OUTPUTFILE}

