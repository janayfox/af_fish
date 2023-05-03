#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#######################################################
### Goal: Run Trinity on EA samples
### Author: Janay Fox
#######################################################

module load singularity/3.8

# singularity exec -e -H /tmp -B `pwd` -B /home/janayfox/scratch/afFishRNA/readsBeforeRmoverrep/BA:/data \
# trinityrnaseq.v2.15.0.simg Trinity --seqType fq --CPU 4 --max_memory 45G --SS_lib_type RF \
# --output /data/trinity_output \
# --left /data/0590g_left.fq,/data/0591g_left.fq,\
# /data/0592g_left.fq,/data/0593g_left.fq,\
# /data/0594g_left.fq,/data/0595g_left.fq,\
# /data/0596g_left.fq,/data/0597g_left.fq,\
# /data/0599g_left.fq,/data/0600g_left.fq,\
# /data/0601g_left.fq,/data/0602g_left.fq,\
# /data/0603g_left.fq,/data/0604g_left.fq,\
# /data/0605g_left.fq,/data/0606g_left.fq,\
# /data/0607g_left.fq,/data/0608g_left.fq,\
# /data/0609g_left.fq,/data/0625g_left.fq,\
# /data/0629g_left.fq,/data/0631g_left.fq,\
# /data/0632g_left.fq,/data/0633g_left.fq,\
# /data/0634g_left.fq,/data/0635g_left.fq,\
# /data/0636g_left.fq,/data/0638g_left.fq,\
# /data/0639g_left.fq,/data/0640g_left.fq,\
# /data/0641g_left.fq,/data/0642g_left.fq,\
# /data/0643g_left.fq,/data/0648g_left.fq,\
# /data/0650g_left.fq,/data/0651g_left.fq,\
# /data/0655g_left.fq \
# --right /data/0590g_right.fq,/data/0591g_right.fq,\
# /data/0592g_right.fq,/data/0593g_right.fq,\
# /data/0594g_right.fq,/data/0595g_right.fq,\
# /data/0596g_right.fq,/data/0597g_right.fq,\
# /data/0599g_right.fq,/data/0600g_right.fq,\
# /data/0601g_right.fq,/data/0602g_right.fq,\
# /data/0603g_right.fq,/data/0604g_right.fq,\
# /data/0605g_right.fq,/data/0606g_right.fq,\
# /data/0607g_right.fq,/data/0608g_right.fq,\
# /data/0609g_right.fq,/data/0625g_right.fq,\
# /data/0629g_right.fq,/data/0631g_right.fq,\
# /data/0632g_right.fq,/data/0633g_right.fq,\
# /data/0634g_right.fq,/data/0635g_right.fq,\
# /data/0636g_right.fq,/data/0638g_right.fq,\
# /data/0639g_right.fq,/data/0640g_right.fq,\
# /data/0641g_right.fq,/data/0642g_right.fq,\
# /data/0643g_right.fq,/data/0648g_right.fq,\
# /data/0650g_right.fq,/data/0651g_right.fq,\
# /data/0655g_right.fq 

singularity exec -e -H /tmp -B `pwd` \
trinityrnaseq.v2.15.0.simg Trinity --seqType fq --CPU 4 --max_memory 45G --SS_lib_type RF \
--output /lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/trinity_output --full_cleanup \
--left /lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0590g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0591g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0592g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0593g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0594g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0595g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0596g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0597g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0599g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0600g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0601g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0602g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0603g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0604g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0605g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0606g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0607g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0608g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0609g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0625g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0629g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0631g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0632g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0633g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0634g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0635g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0636g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0638g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0639g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0640g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0641g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0642g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0643g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0648g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0650g_left.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0651g_left.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0655g_left.fq \
--right /lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0590g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0591g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0592g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0593g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0594g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0595g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0596g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0597g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0599g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0600g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0601g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0602g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0603g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0604g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0605g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0606g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0607g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0608g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0609g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0625g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0629g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0631g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0632g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0633g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0634g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0635g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0636g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0638g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0639g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0640g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0641g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0642g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0643g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0648g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0650g_right.fq,/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0651g_right.fq,\
/lustre04/scratch/janayfox/afFishRNA/readsBeforeRmoverrep/BA/0655g_right.fq 

