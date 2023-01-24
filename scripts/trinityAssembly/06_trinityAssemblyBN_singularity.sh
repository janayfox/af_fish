#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load singularity/3.8

singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BN:/data \
trinityrnaseq.v2.15.0.simg Trinity --seqType fq --CPU 4 --max_memory 45G --SS_lib_type RF \
--output /data/trinity_output \
--left /data/0580g_left.fq,/data/0582g_left.fq,\
/data/0583g_left.fq,/data/0584g_left.fq,\
/data/0585g_left.fq,/data/0586g_left.fq,\
/data/0587g_left.fq,/data/0588g_left.fq,\
/data/0589g_left.fq,/data/0610g_left.fq,\
/data/0611g_left.fq,/data/0612g_left.fq,\
/data/0613g_left.fq,/data/0614g_left.fq,\
/data/0615g_left.fq,/data/0616g_left.fq,\
/data/0617g_left.fq,/data/0619g_left.fq,\
/data/0620g_left.fq,/data/0622g_left.fq,\
/data/0623g_left.fq,/data/0624g_left.fq,\
/data/0626g_left.fq,/data/0627g_left.fq,\
/data/0628g_left.fq,/data/0630g_left.fq,\
/data/0637g_left.fq,/data/0644g_left.fq,\
/data/0645g_left.fq,/data/0646g_left.fq,\
/data/0647g_left.fq,/data/0649g_left.fq,\
/data/0653g_left.fq \
--right /data/0580g_right.fq,/data/0582g_right.fq,\
/data/0583g_right.fq,/data/0584g_right.fq,\
/data/0585g_right.fq,/data/0586g_right.fq,\
/data/0587g_right.fq,/data/0588g_right.fq,\
/data/0589g_right.fq,/data/0610g_right.fq,\
/data/0611g_right.fq,/data/0612g_right.fq,\
/data/0613g_right.fq,/data/0614g_right.fq,\
/data/0615g_right.fq,/data/0616g_right.fq,\
/data/0617g_right.fq,/data/0619g_right.fq,\
/data/0620g_right.fq,/data/0622g_right.fq,\
/data/0623g_right.fq,/data/0624g_right.fq,\
/data/0626g_right.fq,/data/0627g_right.fq,\
/data/0628g_right.fq,/data/0630g_right.fq,\
/data/0637g_right.fq,/data/0644g_right.fq,\
/data/0645g_right.fq,/data/0646g_right.fq,\
/data/0647g_right.fq,/data/0649g_right.fq,\
/data/0653g_right.fq 