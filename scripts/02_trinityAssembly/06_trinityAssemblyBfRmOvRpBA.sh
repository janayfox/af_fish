#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load singularity/3.8

singularity exec -e -B /home/janayfox/scratch/afFishRNA/readsBeforeRmoverrep/BA:/data \
trinityrnaseq.v2.15.0.simg Trinity --seqType fq --CPU 4 --max_memory 45G --SS_lib_type RF \
--output /data/trinity_output \
--left /data/0590g.fq.1,/data/0591g.fq.1,\
/data/0592g.fq.1,/data/0593g.fq.1,\
/data/0594g.fq.1,/data/0595g.fq.1,\
/data/0596g.fq.1,/data/0597g.fq.1,\
/data/0599g.fq.1,/data/0600g.fq.1,\
/data/0601g.fq.1,/data/0602g.fq.1,\
/data/0603g.fq.1,/data/0604g.fq.1,\
/data/0605g.fq.1,/data/0606g.fq.1,\
/data/0607g.fq.1,/data/0608g.fq.1,\
/data/0609g.fq.1,/data/0625g.fq.1,\
/data/0629g.fq.1,/data/0631g.fq.1,\
/data/0632g.fq.1,/data/0633g.fq.1,\
/data/0634g.fq.1,/data/0635g.fq.1,\
/data/0636g.fq.1,/data/0638g.fq.1,\
/data/0639g.fq.1,/data/0640g.fq.1,\
/data/0641g.fq.1,/data/0642g.fq.1,\
/data/0643g.fq.1,/data/0648g.fq.1,\
/data/0650g.fq.1,/data/0651g.fq.1,\
/data/0655g.fq.1 \
--right /data/0590g.fq.2,/data/0591g.fq.2,\
/data/0592g.fq.2,/data/0593g.fq.2,\
/data/0594g.fq.2,/data/0595g.fq.2,\
/data/0596g.fq.2,/data/0597g.fq.2,\
/data/0599g.fq.2,/data/0600g.fq.2,\
/data/0601g.fq.2,/data/0602g.fq.2,\
/data/0603g.fq.2,/data/0604g.fq.2,\
/data/0605g.fq.2,/data/0606g.fq.2,\
/data/0607g.fq.2,/data/0608g.fq.2,\
/data/0609g.fq.2,/data/0625g.fq.2,\
/data/0629g.fq.2,/data/0631g.fq.2,\
/data/0632g.fq.2,/data/0633g.fq.2,\
/data/0634g.fq.2,/data/0635g.fq.2,\
/data/0636g.fq.2,/data/0638g.fq.2,\
/data/0639g.fq.2,/data/0640g.fq.2,\
/data/0641g.fq.2,/data/0642g.fq.2,\
/data/0643g.fq.2,/data/0648g.fq.2,\
/data/0650g.fq.2,/data/0651g.fq.2,\
/data/0655g.fq.2 