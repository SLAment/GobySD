#!/bin/bash
#To use on local machine


#align ddRAD reads from P. flavescens to the new P. falvescens ref genome. 

cd ~/sexdetdir

#Alignment to haplotype 1
mkdir Aligned_hap1
mkdir Aligned_hap1/sorted_hap1
mkdir Aligned_hap1/sorted_hap1/stats

#it takes about 24h on local machine for all gobies (180 individuals)
for var in $(seq -f "%02g" 1 180)
do
~/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 12 ~/sexdetdir/refsequence/refindex/hap1_index ~/realdir/demult_final/final_merged2/DYN-$var".1".fq.gz ~/realdir/demult_final/final_merged2/DYN-$var".2".fq.gz | samtools view --threads 12 -S -b -h - > ./Aligned_hap1/ALH1-$var.bam
samtools sort ./Aligned_hap1/ALH1-$var.bam -o ./Aligned_hap1/sorted_hap1/SORTH1-$var.bam    #sort in order of ref sequence
samtools index ./Aligned_hap1/sorted_hap1/SORTH1-$var.bam	
samtools stats --threads 12 ./Aligned_hap1/sorted_hap1/SORTH1-$var.bam > ./Aligned_hap1/sorted_hap1/stats/SORTH1-$var.stats #create stats file for multiQC
done

multiqc ./Aligned_hap1/sorted_hap1/stats

#Alignment to haplotype 2
mkdir Aligned_hap2
mkdir Aligned_hap2/sorted_hap2
mkdir Aligned_hap2/sorted_hap2/stats

#it takes about 24h on local machine for all gobies
for var in $(seq -f "%02g" 1 180)
do
~/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 12 ~/sexdetdir/refsequence/refindex/hap2_index ~/realdir/demult_final/final_merged2/DYN-$var".1".fq.gz ~/realdir/demult_final/final_merged2/DYN-$var".2".fq.gz | samtools view --threads 12 -S -b -h - > ./Aligned_hap2/ALH2-$var.bam
samtools sort ./Aligned_hap2/ALH2-$var.bam -o ./Aligned_hap2/sorted_hap2/SORTH2-$var.bam    #sort in order of ref sequence
samtools index ./Aligned_hap2/sorted_hap2/SORTH2-$var.bam	
samtools stats --threads 12 ./Aligned_hap2/sorted_hap2/SORTH2-$var.bam > ./Aligned_hap2/sorted_hap2/stats/SORTH2-$var.stats #create stats file for multiQC
done

multiqc ./Aligned_hap2/sorted_hap2/stats
