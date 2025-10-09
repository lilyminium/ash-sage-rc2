#!/bin/bash

host=$(sed 1q host)
port=$(awk '/port/ {print $NF}' optimize.in)

export SLURM_TMPDIR=/tmp
export MYTMPDIR=/tmp/$USER
export TMPDIR=$SLURM_TMPDIR/$SLURM_JOB_NAME

mkdir $TMPDIR
mkdir $MYTMPDIR

worker_num=$(squeue -u $USER | grep wq -c)
ncpus=10

echo submitting worker $worker_num with $ncpus cpus on $host:$port
conda_env=$(sed 1q env.path)

cmd=$(mktemp)
cat << EOF > $cmd
#!/usr/bin/env bash
#SBATCH -J wq-$port
#SBATCH -p free
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=${ncpus}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --array=1-300%30
#SBATCH --account lilyw7
# SBATCH --export ALL
#SBATCH -o worker-logs/worker-${worker_num}-%a.log

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

source $HOME/.bashrc

mkdir ${MYTMPDIR} -p
cd $MYTMPDIR

if [[ ! -d $conda_env ]]; then
    compressed_env=/dfs9/dmobley-lab/$USER/envs/$conda_env.tar.gz
    cp $compressed_env .
    mkdir -p $conda_env
    tar xzf $compressed_env -C $conda_env
fi

conda activate $conda_env

sleep 300

for i in \$(seq  \$SLURM_NTASKS ); do
        echo $i
        work_queue_worker --cores 1 -s ${MYTMPDIR} \
                          --disk-threshold=0.002 --disk=3000 \
                          --memory-threshold=1000 -t 3600 -b 20 \
                          --memory=1000 $host:$port &
done
wait
EOF

# wait 10 min

sbatch $@ $cmd
rm $cmd
