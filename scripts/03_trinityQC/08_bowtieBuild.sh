#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

#$1 = left handed reads
#$2 = right haded reads
#$3 = sample name

module load nixpkgs/16.09
module load gcc/7.3.0
module load intel/2018.3
module load bowtie2/2.3.4.3
#module load StdEnv/2020
#module load samtools/1.16.1

#build index BN
bowtie2-build ./BA.Trinity.fasta BA.Trinity.fasta

#build index BA
bowtie2-build ./BN.Trinity.fasta BN.Trinity.fasta


#run alignment 
bowtie2 -p 10 -q --no-unal -k 20 -x BA.Trinity.fasta -1 $1 -2 $2 2>align_stats.txt
#samtools view -@10 -Sb -o ${3}_bowtie2.bam

#for BA
sbatch 07_bowtieBuild.sh \
0590g_left.fq,0591g_left.fq,\
0592g_left.fq,0593g_left.fq,\
0594g_left.fq,0595g_left.fq,\
0596g_left.fq,0597g_left.fq,\
0599g_left.fq,0600g_left.fq,\
0601g_left.fq,0602g_left.fq,\
0603g_left.fq,0604g_left.fq,\
0605g_left.fq,0606g_left.fq,\
0607g_left.fq,0608g_left.fq,\
0609g_left.fq,0625g_left.fq,\
0629g_left.fq,0631g_left.fq,\
0632g_left.fq,0633g_left.fq,\
0634g_left.fq,0635g_left.fq,\
0636g_left.fq,0638g_left.fq,\
0639g_left.fq,0640g_left.fq,\
0641g_left.fq,0642g_left.fq,\
0643g_left.fq,0648g_left.fq,\
0650g_left.fq,0651g_left.fq,\
0655g_left.fq \
0590g_right.fq,0591g_right.fq,\
0592g_right.fq,0593g_right.fq,\
0594g_right.fq,0595g_right.fq,\
0596g_right.fq,0597g_right.fq,\
0599g_right.fq,0600g_right.fq,\
0601g_right.fq,0602g_right.fq,\
0603g_right.fq,0604g_right.fq,\
0605g_right.fq,0606g_right.fq,\
0607g_right.fq,0608g_right.fq,\
0609g_right.fq,0625g_right.fq,\
0629g_right.fq,0631g_right.fq,\
0632g_right.fq,0633g_right.fq,\
0634g_right.fq,0635g_right.fq,\
0636g_right.fq,0638g_right.fq,\
0639g_right.fq,0640g_right.fq,\
0641g_right.fq,0642g_right.fq,\
0643g_right.fq,0648g_right.fq,\
0650g_right.fq,0651g_right.fq,\
0655g_right.fq \
0590g,0591g,\
0592g,0593g,\
0594g,0595g,\
0596g,0597g,\
0599g,0600g,\
0601g,0602g,\
0603g,0604g,\
0605g,0606g,\
0607g,0608g,\
0609g,0625g,\
0629g,0631g,\
0632g,0633g,\
0634g,0635g,\
0636g,0638g,\
0639g,0640g,\
0641g,0642g,\
0643g,0648g,\
0650g,0651g,\
0655g

# for BN 
sbatch 07_bowtieBuild.sh \
0580g_left.fq,0582g_left.fq,\
0583g_left.fq,0584g_left.fq,\
0585g_left.fq,0586g_left.fq,\
0587g_left.fq,0588g_left.fq,\
0589g_left.fq,0610g_left.fq,\
0611g_left.fq,0612g_left.fq,\
0613g_left.fq,0614g_left.fq,\
0615g_left.fq,0616g_left.fq,\
0617g_left.fq,0619g_left.fq,\
0620g_left.fq,0622g_left.fq,\
0623g_left.fq,0624g_left.fq,\
0626g_left.fq,0627g_left.fq,\
0628g_left.fq,0630g_left.fq,\
0637g_left.fq,0644g_left.fq,\
0645g_left.fq,0646g_left.fq,\
0647g_left.fq,0649g_left.fq,\
0653g_left.fq \
0580g_right.fq,0582g_right.fq,\
0583g_right.fq,0584g_right.fq,\
0585g_right.fq,0586g_right.fq,\
0587g_right.fq,0588g_right.fq,\
0589g_right.fq,0610g_right.fq,\
0611g_right.fq,0612g_right.fq,\
0613g_right.fq,0614g_right.fq,\
0615g_right.fq,0616g_right.fq,\
0617g_right.fq,0619g_right.fq,\
0620g_right.fq,0622g_right.fq,\
0623g_right.fq,0624g_right.fq,\
0626g_right.fq,0627g_right.fq,\
0628g_right.fq,0630g_right.fq,\
0637g_right.fq,0644g_right.fq,\
0645g_right.fq,0646g_right.fq,\
0647g_right.fq,0649g_right.fq,\
0653g_right.fq \
0580g,0582g,\
0583g,0584g,\
0585g,0586g,\
0587g,0588g,\
0589g,0610g,\
0611g,0612g,\
0613g,0614g,\
0615g,0616g,\
0617g,0619g,\
0620g,0622g,\
0623g,0624g,\
0626g,0627g,\
0628g,0630g,\
0637g,0644g,\
0645g,0646g,\
0647g,0649g,\
0653g 
