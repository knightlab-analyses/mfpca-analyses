#!/bin/bash

#PBS -N mfpca_t2d_new
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=4
#PBS -l pmem=10gb
#PBS -V
#PBS -j oe
#PBS -d .
#PBS -t 3
#PBS -o ../messages_outputs/
#PBS -e ../messages_errors/

set -e
cpus=$PBS_NUM_PPN

export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER/mfpca-analyses
[ ! -d $TMPDIR ] && mkdir $TMPDIR
export TMPDIR=$TMPDIR/applications/HMP_t2d

# do something
module load R_3.5.3-9.3.0 binutils_2.30

# t: number of candidate models (e.g. 1-7%10)
# fit MFPCA all on candidate models
#Rscript 1_mfpca_t2d_microb_protein_cytok_a1c_model_selection.R $TMPDIR

# t: total number of candidate models (e.g. 7)
# compare LOOIC on all models
Rscript 2_mfpca_t2d_model_optimal.R $TMPDIR
