#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=1-00:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load StdEnv/2020 
module load gcc/9.3.0
module load blast+/2.13.0
module load singularity/3.8

#build a blastable database 
makeblastdb -in uniprot_sprot.fasta -dbtype prot

#perform blast search, reporting only the top alignment:
blastx -query BA.Trinity.fasta -db uniprot_sprot.fasta -out BA.blastx.outfmt6 \
-evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6

blastx -query BN.Trinity.fasta -db uniprot_sprot.fasta -out BN.blastx.outfmt6 \
-evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6

#examine the percentage of the target being aligned to the best matching Trinity transcript
singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BA:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/analyze_blastPlus_topHit_coverage.pl \
BA.blastx.outfmt6 /data/BA.Trinity.fasta /data/uniprot_sprot.fasta

singularity exec -e -B /home/janayfox/scratch/afFishRNA/cleanedReads/BN:/data \
trinityrnaseq.v2.15.0.simg /usr/local/bin/util/analyze_blastPlus_topHit_coverage.pl \
BN.blastx.outfmt6 /data/BN.Trinity.fasta /data/uniprot_sprot.fasta