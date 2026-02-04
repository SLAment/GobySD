#!/bin/bash
#To use on local machine 


cd ~/sexdetdir
mkdir catalogg_hap1
mkdir catalogg_hap2

gstacks -t 12 -I ./Aligned_hap1/sorted_hap1 -M ./popmap_hap1.txt -O ./catalogg_hap1  #bamfiles have to be sorted or error message
gstacks -t 12 -I ./Aligned_hap2/sorted_hap2 -M ./popmap_hap2.txt -O ./catalogg_hap2  #bamfiles have to be sorted or error message