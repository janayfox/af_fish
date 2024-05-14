#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-01:00
#SBATCH --account=def-barrett
#SBATCH --mail-type=ALL
#SBATCH --mail-user=janay.fox@mail.mcgill.ca

module load StdEnv/2020
module load gcc/9.3.0
#module load StdEnv/2023
#module load nixpkgs/16.09 
#module load intel/2018.3
#module load gatk/4.4.0.0 
#module load plink/2.00 
#module load vcftools/0.1.16 
module load bcftools/1.16
#module load picard/3.1.0

## new version ## 
#merge VCFs
# java -jar $EBROOTPICARD/picard.jar MergeVcfs \
# I=BA.list.vcf.txt \
# O=BA_variants.vcf 

# java -jar $EBROOTPICARD/picard.jar MergeVcfs \
# I=BN.list.vcf.txt \
# O=BN_variants.vcf  

#merge using bcftools 
# cat BA.list.vcf.txt| while read file 
# do
# bgzip -c "$file" > "${file}.gz"
# done

# Sort and index VCF files listed in BA.list.vcf.txt
while read -r file; do
    # Sort VCF file
    sorted_file="${file%.vcf}.sorted.vcf.gz"
    bcftools sort "$file" -Oz -o "$sorted_file"

    # Index sorted VCF file
    bcftools index "$sorted_file"
done < BA.list.vcf.txt

# Sort and index VCF files listed in BN.list.vcf.txt
while read -r file; do
    # Sort VCF file
    sorted_file="${file%.vcf}.sorted.vcf.gz"
    bcftools sort "$file" -Oz -o "$sorted_file"

    # Index sorted VCF file
    bcftools index "$sorted_file"
done < BN.list.vcf.txt

#bcftools merge –file-list BA.list.vcf.txt -Oz -o BA_merged_vcf_filt.vcf.gz
#bcftools merge –file-list BN.list.vcf.txt -Oz -o BN_merged_vcf_filt.vcf.gz

# #filter sites with missing data, minor allele frequency < 5% and out of hwe, 
#vcftools --vcf BA_variants.vcf --max-alleles 2 --max-missing 1.0 --maf 0.05 --hwe 0.05 --recode --recode-INFO-all --out BA_variants_filtered.vcf
#vcftools --vcf BN_variants.vcf --max-alleles 2 --max-missing 1.0 --maf 0.05 --hwe 0.05 --recode --recode-INFO-all --out BN_variants_filtered.vcf

# #filter linkage disequilibirum
# bcftools +prune -m 0.8 -w 1000 BA_variants_filtered.vcf -Ov -o BA_variants_filtered_ld.vcf
# bcftools +prune -m 0.8 -w 1000 BN_variants_filtered.vcf -Ov -o BN_variants_filtered_ld.vcf
