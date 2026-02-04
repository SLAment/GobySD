#!/bin/bash
#take vcf from Stacks pipeline after alignment to our P. flavescens genome(s) and extract sex associated SNPs by frequency and missingness

cd ~/sexdetdir


#vcf filtration stage
#vcftools --vcf outpop_hap1/populations.snps.vcf --min-alleles 2 --max-alleles 2 --min-meanDP 5  --max-meanDP 40 --recode --out sexassoc/vcf_hap1
#vcftools --vcf outpop_hap2/populations.snps.vcf --min-alleles 2 --max-alleles 2 --min-meanDP 5  --max-meanDP 40 --recode --out sexassoc/vcf_hap2

cd sexassoc
#make bed file for plink
#~/plink_files/plink --vcf vcf_hap1.recode.vcf --make-bed --allow-extra-chr --out plink_hap1
#~/plink_files/plink --vcf vcf_hap2.recode.vcf --make-bed --allow-extra-chr --out plink_hap2

#1.
#test association of allele frequency with sex
#~/plink_files/plink --bfile plink_hap1 --pheno phenotypes.txt --pheno-name SEX --assoc fisher --allow-extra-chr --allow-no-sex --out assocsex_hap1
~/plink_files/plink --bfile plink_hap2 --pheno phenotypes2.txt --pheno-name SEX --assoc fisher --allow-extra-chr --allow-no-sex --out assocsex_hap2
#switch to R for visual inspection of frequencies and selection of candidat SNPs.

#2.
#test association of locus missingeness with sex
#~/plink_files/plink --bfile plink_hap1 --pheno phenotypes.txt --pheno-name SEX --test-missing --allow-extra-chr --allow-no-sex --out missingsex_hap1
~/plink_files/plink --bfile plink_hap2 --pheno phenotypes2.txt --pheno-name SEX --test-missing --allow-extra-chr --allow-no-sex --out missingsex_hap2