#!/bin/bash

#PBS -N mfpca_t2d_new
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=4
#PBS -l pmem=10gb
#PBS -V
#PBS -j oe
#PBS -d .
#PBS -t 1
#PBS -o ../messages_outputs/
#PBS -e ../messages_errors/

set -e
cpus=$PBS_NUM_PPN

export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER/Oct_2020_mfpca
[ ! -d $TMPDIR ] && mkdir $TMPDIR
export TMPDIR=$TMPDIR/applications/HMP_t2d

# do something
module load R_3.5.3-9.3.0 binutils_2.30

# t: number of candidate models (e.g. 1-7%10)
# fit MFPCA all on candidate models
#Rscript 1_mfpca_t2d_microb_protein_cytok_model_selection.R $TMPDIR
#Rscript 2_1_mfpca_t2d_A1c_PDGFBB_IGHA2_ratio_v2.R $TMPDIR


# t: total number of candidate models (e.g. 7)
# compare LOOIC on all models
# Rscript 1_2_mfpca_t2d_model_optimal.R $TMPDIR
#Rscript 2_2_mfpca_t2d_model_optimal_v2.R $TMPDIR


# t: optimal model index (e.g. 1)
# model diagnostics + results visualization on optimal model
#Rscript 1_3_mfpca_t2d_model_results.R $TMPDIR


# compare computational time
#Rscript 3_1_computation_time_mod3.R $TMPDIR
#Rscript 3_2_computation_time_mod5.R $TMPDIR
Rscript 3_3_computation_time_mod7.R $TMPDIR
