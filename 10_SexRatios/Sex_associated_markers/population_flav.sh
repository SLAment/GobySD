#!/bin/bash
#To use on local machine 


cd ~/sexdetdir
mkdir outpop_hap1
mkdir outpop_hap2


# -p minimum number of populations a locus must be present in to be kept . 1= keep all
# -r minium percentage of individuals per pop to have the locus to keep it.
# -R minimum percentage accross pops. These should all disappear when I filter for missing data later anyways... I could be more stringent

populations -t 12 --structure --vcf --fstats -P ./catalogg_hap1 -M ./popmap_hap1.txt -O ./outpop_hap1 -p 1 -r 0.1 -R 0.1

populations -t 12 --structure --vcf --fstats -P ./catalogg_hap2 -M ./popmap_hap2.txt -O ./outpop_hap2 -p 1 -r 0.1 -R 0.1