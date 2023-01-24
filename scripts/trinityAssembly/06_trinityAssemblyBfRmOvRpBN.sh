#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=45G
#SBATCH --time=7-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load singularity/3.8

singularity exec -e -B /home/janayfox/scratch/afFishRNA/readsBeforeRmoverrep/BN:/data \
trinityrnaseq.v2.15.0.simg Trinity --seqType fq --CPU 4 --max_memory 45G --SS_lib_type RF \
--output /data/trinity_output \
--left /data/0580g.fq.1,/data/0582g.fq.1,\
/data/0583g.fq.1,/data/0584g.fq.1,\
/data/0585g.fq.1,/data/0586g.fq.1,\
/data/0587g.fq.1,/data/0588g.fq.1,\
/data/0589g.fq.1,/data/0610g.fq.1,\
/data/0611g.fq.1,/data/0612g.fq.1,\
/data/0613g.fq.1,/data/0614g.fq.1,\
/data/0615g.fq.1,/data/0616g.fq.1,\
/data/0617g.fq.1,/data/0619g.fq.1,\
/data/0620g.fq.1,/data/0622g.fq.1,\
/data/0623g.fq.1,/data/0624g.fq.1,\
/data/0626g.fq.1,/data/0627g.fq.1,\
/data/0628g.fq.1,/data/0630g.fq.1,\
/data/0637g.fq.1,/data/0644g.fq.1,\
/data/0645g.fq.1,/data/0646g.fq.1,\
/data/0647g.fq.1,/data/0649g.fq.1,\
/data/0653g.fq.1 \
--right /data/0580g.fq.2,/data/0582g.fq.2,\
/data/0583g.fq.2,/data/0584g.fq.2,\
/data/0585g.fq.2,/data/0586g.fq.2,\
/data/0587g.fq.2,/data/0588g.fq.2,\
/data/0589g.fq.2,/data/0610g.fq.2,\
/data/0611g.fq.2,/data/0612g.fq.2,\
/data/0613g.fq.2,/data/0614g.fq.2,\
/data/0615g.fq.2,/data/0616g.fq.2,\
/data/0617g.fq.2,/data/0619g.fq.2,\
/data/0620g.fq.2,/data/0622g.fq.2,\
/data/0623g.fq.2,/data/0624g.fq.2,\
/data/0626g.fq.2,/data/0627g.fq.2,\
/data/0628g.fq.2,/data/0630g.fq.2,\
/data/0637g.fq.2,/data/0644g.fq.2,\
/data/0645g.fq.2,/data/0646g.fq.2,\
/data/0647g.fq.2,/data/0649g.fq.2,\
/data/0653g.fq.2 
