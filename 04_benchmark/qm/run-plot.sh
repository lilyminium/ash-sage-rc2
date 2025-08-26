#!/usr/bin/env bash
#SBATCH -J plot
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

export OE_LICENSE="/data/homezvol3/lilyw7/oe_license.txt"

mkdir images 

python plot-ddes.py \
    -ff 'Parsley 1.3.1' 'openff_unconstrained-1.3.1' \
    -ff 'Sage 2.0.0' "openff_unconstrained-2.0.0" \
    -ff 'Sage 2.2.1' "openff_unconstrained-2.2.1" \
    -ff 'Sage 2.2.1 + AshGC' "openff_unconstrained-2.2.1-ashgc" \
    -ff 'Sage 2.3.0rc1' "openff_unconstrained-2.3.0rc1" \
    -ff 'v1-k100' "fb-fit-v1-single-mean-k100_unconstrained" \
    -ff 'v3-k100' "fb-fit-v3-single-mean-k100_unconstrained" \
    -i  ddes  -o  images


python plot-rmsd-tfd.py \
    -ff 'Parsley 1.3.1' 'openff_unconstrained-1.3.1' \
    -ff 'Sage 2.0.0' "openff_unconstrained-2.0.0" \
    -ff 'Sage 2.2.1' "openff_unconstrained-2.2.1" \
    -ff 'Sage 2.2.1 + AshGC' "openff_unconstrained-2.2.1-ashgc" \
    -ff 'Sage 2.3.0rc1' "openff_unconstrained-2.3.0rc1" \
    -ff 'v1-k100' "fb-fit-v1-single-mean-k100_unconstrained" \
    -ff 'v3-k100' "fb-fit-v3-single-mean-k100_unconstrained" \
    -i rmsd-tfd -o images
