#!/bin/bash
#SBATCH -J fit-ff
#SBATCH -p standard
#SBATCH -t 144:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=10000mb
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mail-user=lilyw7@uci.edu
#SBATCH --constraint=fastscratch
#SBATCH -o master-%A.out
#SBATCH -e master-%A.err

source $HOME/.bashrc
conda_env=nagl-valence
conda activate $conda_env

# set up worker run
echo $(hostname) > $SLURM_SUBMIT_DIR/host
echo $CONDA_PREFIX > $SLURM_SUBMIT_DIR/env.path
bash submit_hpc3_worker_local.sh

# set up and move to tmpdir
TMPDIR=/tmp/$USER/$SLURM_JOB_ID
rm -rf $TMPDIR
mkdir -p $TMPDIR
cd $TMPDIR

# move conda environment over
compressed_env=/dfs9/dmobley-lab/$USER/envs/$conda_env.tar.gz
cp $compressed_env .
mkdir -p $conda_env
tar xzf $compressed_env -C $conda_env

# copy files over
scp -C  $SLURM_SUBMIT_DIR/optimize.in     $TMPDIR
scp -C  $SLURM_SUBMIT_DIR/targets.tar.gz  $TMPDIR
scp -Cr $SLURM_SUBMIT_DIR/forcefield      $TMPDIR

tar -xzf targets.tar.gz

datadir=$(pwd)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1


# in case of ownership issues messing with disk quota
newgrp dmobley_lab_share
mkdir -p result/optimize
touch result/optimize/force-field.offxml

if ForceBalance.py optimize.in ; then
   newgrp dmobley_lab_share
   echo "-- Force field done --"
   echo "-- -- -- -- -- -- -- --"
   tar -czf optimize.tmp.tar.gz optimize.tmp
   rsync  -azIi -rv --exclude="optimize.tmp" --exclude="optimize.bak" \
	  --exclude="fb*" --exclude="nagl*"\
	  --chown=lilyw7:dmobley_lab_share  \
	  --exclude="targets*" $TMPDIR/* $SLURM_SUBMIT_DIR > copy.log
   rm -rf $TMPDIR
fi

echo "All done"

cd $SLURM_SUBMIT_DIR

# cancel rest of CPU jobs
PORT=$(awk '/port/ {print $NF}' optimize.in)
WORKER_ID=$(squeue -u lilyw7 | grep $PORT | awk -F'_' '{print $1}' | head -n 1)
echo "cancelling $WORKER_ID"
scancel $WORKER_ID


tar -czvf worker-logs.tar.gz worker-logs
rm -rf worker-logs
