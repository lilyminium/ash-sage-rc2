#!/usr/bin/env bash
#SBATCH -J sdf
#SBATCH -p standard
#SBATCH -t 1:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate nagl-valence

# QCA_ID=36973260
QCA_ID=36984578
QCA_ID=138534827
QCA_ID=120165274
QCA_ID=120192665
QCA_ID=137157615
QCA_ID=120191583
QCA_ID=36966336
QCA_ID=146498204
QCA_ID=120167857
QCA_ID=36973260
QCA_ID=36973258
QCA_ID=36973257
QCA_ID=6098249
QCA_ID=120191135
QCA_ID=137149024 # bonds
QCA_ID=36990645 # torsions
QCA_ID=3065532 # torsions 180
QCA_ID=138556169 # torsions 180
QCA_ID=36957427 #b77

# t145, t157, t148, t149 conformer geometries
QCA_ID=137149023
QCA_ID=137149024
QCA_ID=137149025
QCA_ID=137149026
QCA_ID=137149027
QCA_ID=137149028

#t151a
QCA_ID=137149093

mkdir -p sdfs

python retrieve-sdf-for-qca-id.py   \
    -i  data/optimization           \
    -o  sdfs/qca-id-${QCA_ID}.sdf   \
    -q  ${QCA_ID}

