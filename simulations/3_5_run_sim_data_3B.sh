#!/bin/bash

#PBS -N sim_mfpca
#PBS -l walltime=40:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=10gb
#PBS -V
#PBS -j oe
#PBS -d .
#PBS -t 500
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
#Rscript 3_0_sim_3B.R $TMPDIR # t = 1 # export TMPDIR=$TMPDIR/simulations/sim_data_NB (need to correct output directory to match other code

#Rscript 3_1_estimate_all_3B.R $TMPDIR # t = 1-1000%10
#Rscript 3_1_estimate_all_3B_sub.R $TMPDIR # t = 11-1000%10 (check what's needed)
#Rscript 3_1_estimate_all_3B_sub2.R $TMPDIR # t = 1-100%10


#Rscript 3_2_estimate_all_3B_process.R $TMPDIR # t = 1-1000%10
Rscript 3_3_estimate_all_3B_summary.R $TMPDIR # t = 1000 (test t = 500)
#Rscript 3_4_estimate_all_computation_time.R $TMPDIR # t = 1000

#Rscript 4_0_sim_5B.R $TMPDIR # t = 1 # export TMPDIR=$TMPDIR/simulations
#Rscript 4_1_estimate_5B.R $TMPDIR # t = 1-1000%10

#Rscript 5_0_sim_10B.R $TMPDIR # t = 1
#Rscript 5_1_estimate_10B.R $TMPDIR # t = 1-100%10

#mv $tmp/outdir ./outdir

