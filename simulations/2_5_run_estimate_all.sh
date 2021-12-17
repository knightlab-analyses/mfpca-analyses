#!/bin/bash

#PBS -N sim_mfpca
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=10gb
#PBS -V
#PBS -j oe
#PBS -d .
#PBS -t 1000
#PBS -o messages_outputs/
#PBS -e messages_errors/

set -e
cpus=$PBS_NUM_PPN

export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER/Oct_2020_mfpca
[ ! -d $TMPDIR ] && mkdir $TMPDIR
export TMPDIR=$TMPDIR/simulations
[ ! -d $TMPDIR ] && mkdir $TMPDIR
#tmp=$(mktemp -d --tmpdir)
#export TMPDIR=$tmp
#trap "rm -r $tmp; unset TMPDIR" EXIT

# do something
module load R_3.5.3-9.3.0 binutils_2.30

# t: number of each simulation (e.g. 1-1000%10)
# estimate all parameters (with default weak prior)
# Rscript 2_0_1_estimate_all_v0_weakPrior.R $TMPDIR
# Rscript 2_1_1_estimate_all_v1_weakPrior.R $TMPDIR
# Rscript 2_2_1_estimate_all_v2_weakPrior.R $TMPDIR
# Rscript 2_3_1_estimate_all_v3_weakPrior.R $TMPDIR
#Rscript 2_4_1_estimate_all_v4_weakPrior.R $TMPDIR

# t: number of each simulation (e.g. 1-1000%10)
#Rscript 2_0_2_estimate_all_v0_weakPrior_process.R $TMPDIR
#Rscript 2_1_2_estimate_all_v1_weakPrior_process.R $TMPDIR
#Rscript 2_2_2_estimate_all_v2_weakPrior_process.R $TMPDIR
#Rscript 2_3_2_estimate_all_v3_weakPrior_process.R $TMPDIR
#Rscript 2_4_2_estimate_all_v4_weakPrior_process.R $TMPDIR

# t: number of total simulations (e.g. 1000)
# Rscript 2_0_3_estimate_all_v0_weakPrior_summary.R $TMPDIR
# Rscript 2_1_3_estimate_all_v1_weakPrior_summary.R $TMPDIR
# Rscript 2_2_3_estimate_all_v2_weakPrior_summary.R $TMPDIR
# Rscript 2_3_3_estimate_all_v3_weakPrior_summary.R $TMPDIR
Rscript 2_4_3_estimate_all_v4_weakPrior_summary.R $TMPDIR

