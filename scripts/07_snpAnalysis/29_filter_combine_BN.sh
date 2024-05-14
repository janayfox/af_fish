#!/bin/bash
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=50G
#SBATCH --time=0-00:05
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load gatk/4.4.0.0 
#module load plink/2.00 
#module load vcftools/0.1.16 
#module load bcftools/1.16

# gatk --java-options "-Xmx300G" CombineGVCFs -R /home/janayfox/scratch/af_fish_RNA/super_transcripts/BA/BA_super_trinity_genes.fasta \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0590g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0591g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0592g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0593g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0594g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0595g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0596g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0597g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0599g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0600g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0601g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0602g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0603g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0604g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0605g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0606g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0607g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0608g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0609g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0625g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0629g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0631g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0632g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0633g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0634g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0635g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0636g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0638g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0639g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0640g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0641g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0642g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0643g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0648g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0650g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0651g/output.vcf \
# -V /home/janayfox/scratch/af_fish_RNA/variant_calling/BA/0655g/output.vcf \
# -O BA.combined.vcf

gatk --java-options "-Xmx150G" CombineGVCFs -R /home/janayfox/scratch/af_fish_RNA/super_transcripts/BN/BN_super_trinity_genes.fasta \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0580g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0582g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0583g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0584g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0585g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0586g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0587g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0588g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0589g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0610g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0611g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0612g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0613g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0614g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0615g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0617g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0619g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0620g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0622g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0623g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0624g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0626g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0627g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0628g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0630g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0637g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0644g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0645g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0646g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0647g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0649g/output.vcf \
-V /home/janayfox/scratch/af_fish_RNA/variant_calling/BN/0653g/output.vcf \
-O BN.combined.vcf

#does genotyping on a single input which is the output from CombineGVCFs
gatk --java-options "-Xmx300G"  GenotypeGVCFs \
-R /home/janayfox/scratch/af_fish_RNA/super_transcripts/BN/BN_super_trinity_genes.fasta \
-O BN.geno.combined.vcf \
-V BN.combined.vcf

# # maybe have to do to get only SNPs? but I dont think so -> gatk SelectVariants

#filter fisher strand and qualbydepth
gatk VariantFiltration \
-R /home/janayfox/scratch/af_fish_RNA/super_transcripts/BN/BN_super_trinity_genes.fasta \
-V BN.geno.combined.vcf \
-O BN.filt.geno.combined.vcf \
--filter-expression "FS > 30.0" --filter-name "FS30" \
--filter-expression "QD < 2.0" --filter-name "QD2" \

#filter sites with missing data, minor allele frequency < 5% and out of hwe, 
vcftools --vcf BN.filt.geno.combined.vcf --max-alleles 2 --max-missing 1 --maf 0.05 --hwe 0.05 --recode --recode-INFO-all --out BN.filt.2.geno.combined.vcf

#filter linkage disequilibirum
bcftools +prune -l 0.8 -w 1000 BN.filt.2.geno.combined.vcf -Ov -o BN.filt.final.geno.combined.vcf