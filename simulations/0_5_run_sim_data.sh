#!/bin/bash

#PBS -N sim_mfpca
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=10gb
#PBS -V
#PBS -j oe
#PBS -d .
#PBS -o messages_outputs/
#PBS -e messages_errors/

set -e
cpus=$PBS_NUM_PPN

export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER/Oct_2020_mfpca
[ ! -d $TMPDIR ] && mkdir $TMPDIR
export TMPDIR=$TMPDIR/simulations/sim_data
[ ! -d $TMPDIR ] && mkdir $TMPDIR
#tmp=$(mktemp -d --tmpdir)
#export TMPDIR=$tmp
#trap "rm -r $tmp; unset TMPDIR" EXIT

# do something
module load R_3.5.3-9.3.0 binutils_2.30
#Rscript 0_0_sim_N100_C80_3B_v0.R $TMPDIR
#Rscript 0_1_sim_N100_C80_3B_v1.R $TMPDIR
#Rscript 0_2_sim_N100_C80_3B_v2.R $TMPDIR
#Rscript 0_3_sim_N100_C80_3B_v3.R $TMPDIR
Rscript 0_4_sim_N100_C80_3B_v4.R $TMPDIR

#mv $tmp/outdir ./outdir

