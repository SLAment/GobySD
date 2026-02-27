#!/bin/bash

mkdir ~/sexdetdir/refsequence/refindex

#Index the desired reference genomes for later alignment with bwamem (time for 1 genome on local machine: 6 minutes)
~/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p hap1_index ~/sexdetdir/refsequence/hap1.p_ctg.fasta
~/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p hap2_index ~/sexdetdir/refsequence/hap2.p_ctg.fasta
#~/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p primary_index ~/sexdetdir/refsequence/p_ctg.fasta