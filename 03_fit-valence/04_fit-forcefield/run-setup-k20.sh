#!/bin/bash
#SBATCH -J generate-fit-k20
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
# conda activate nagl-valence # bespokefit requires pydantic 1
conda activate fb_196_ic_0326

ffver="v0"
ffver="v1"
#ffver="v2"
#ffver="v3"

# agg="median"
agg="mean"

conf="multi"
conf="single"


echo "==== Generating FF for ${ffver} with ${agg} MSM, ${conf} conformers ===="

TAG="fb-fit-${ffver}-${conf}-${agg}-k20"
echo $TAG

FF_FILE="../03_generate-msm/output/msm-force-field-${ffver}-${agg}.offxml"

# chemper currently cannot handle `]~1` type SMIRKS so the FF needs temporary editing
# see issues below
# https://github.com/MobleyLab/chemper/issues/99
# https://github.com/openforcefield/openff-bespokefit/issues/373
# note: not applicable to v0, v1
FF_FILE="msm-force-field-${ffver}-${agg}.offxml"
sed 's/\]~1"/]1"/g' "../03_generate-msm/output/${FF_FILE}" > $FF_FILE


python create-fb-inputs-nagl.py                                                                         \
    --tag                       $TAG                                                                    \
    --optimization-dataset      "../02_curate-data/output/optimizations-${conf}-${ffver}.json"          \
    --torsion-dataset           "../02_curate-data/output/torsiondrives-${ffver}.json"                  \
    --forcefield                $FF_FILE      \
    --valence-counts            "../02_curate-data/counts/valence-counts-${conf}-${ffver}.json"         \
    --torsion-counts            "../02_curate-data/counts/torsion-counts-${ffver}.json"                 \
    --n-min-valence             3   \
    --n-min-torsion             1   \
    --frozen-angle-file         "../03_generate-msm/linear-angles.json"                      \
    --max-iterations            100                                                                 \
    --port                      59980                                                               \
    --output-directory          "output"                                                            \
    --verbose

cd $TAG

cp forcefield/force-field.offxml tmp.offxml

# now need to reverse the bond change
sed 's/\]1"/]~1"/g' tmp.offxml > forcefield/force-field.offxml

tar -czvf targets.tar.gz targets
rm -rf targets

echo "done"

# copy over to different filesystem and run
# this is just a hpc3 quirk
cd ..
cp -r $TAG /pub/lilyw7/ash-sage/
cd /pub/lilyw7/ash-sage/$TAG
cp /dfs9/dmobley-lab/lilyw7/ash-sage/03_fit-valence/04_fit-forcefield/runfiles/* .
sbatch run-fit.sh

pwd

date
