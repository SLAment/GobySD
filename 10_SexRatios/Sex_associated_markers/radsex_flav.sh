#!/bin/bash

#Rad sex, uding RARseq data to indentify sex specific sequences
#Now that we have our reference genome(s), see what we find.
#Radsex only works with forward reads. Start by moving the forward reads to a separate directory
##After moving the forward reads to "forwardreads" I need to remove the ".1" that indicates forward, so that later the population map isn't confused
cd ~/realdir/demult_final/forwardreads
for i in *.1.fq.gz
do mv -- "$i" "${i%.1.fq.gz}.fq.gz"
done

#Process the read files
radsex process -d 1 --threads 12 --input-dir ~/realdir/demult_final/forwardreads --output-file markers_table1_out.tsv
#Find markers associated with sex
radsex signif --markers-table markers_table1_out.tsv --output-file markers1_out.tsv --popmap popmap_radsex.tsv --min-depth 1 --groups m,f
#radsex subset --markers-table realdir/Radsex_files/markers_table1.tsv --output-file realdir/Radsex_files/markers_sextest.tsv --popmap realdir/Radsex_files/popmap_radsex.tsv --groups m,f --min-group1 25 --max-group2 15


#hap1 mapping
~/radsex/radsex-1.2.0/bin/radsex map --markers-file realdir/Radsex_files/markers_table1.tsv --output-file sexdetdir/mapsex_hap1.tsv --popmap realdir/Radsex_files/popmap_radsex.tsv --genome-file sexdetdir/refsequence/hap1.p_ctg.fasta --min-depth 2 --groups m,f
#hap2
~/radsex/radsex-1.2.0/bin/radsex map --markers-file realdir/Radsex_files/markers_table1.tsv --output-file realdir/Radsex_files/mapsex_hap2_DP2.tsv --popmap realdir/Radsex_files/popmap_radsex_out.tsv --genome-file sexdetdir/refsequence/hap2.p_ctg.fasta --min-depth 2 --groups m,f
#~/radsex/radsex-1.2.0/bin/radsex map --markers-file realdir/Radsex_files/markers_table1_out.tsv --output-file realdir/Radsex_files/mapsex_hap2_out.tsv --popmap realdir/Radsex_files/popmap_radsex_out.tsv --genome-file sexdetdir/refsequence/hap2.p_ctg.fasta --min-depth 1 --groups m,f

