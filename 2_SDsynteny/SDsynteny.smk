# -*- snakemake -*-

### SDsynteny
#############################################################################
# Comparing the synteny and annotation of the SD scaffolds in the two-spotted goby *Pomatoschistus flavescens*
# Based on PodoSynteny.smk
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/11/04 
# ++++++++++++++++++++++++++++++++++++++++++++++

from Bio import SeqIO

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"

# Input data
ASSEMBLY = config["assembly"]
satDNAgff = config["satDNAgff"]
SDscaffolds = config["SDscaffolds"]
genegff = config["genegff"]

# Scripts
gggenomes = config["gggenomes"]

# Filtering parameters
MINSIZE = config["MINSIZE"]				# Min length of a hit to survive filtering	        # Min length of a hit to survive filtering for the 3end
EXTRABPSIDE = config["EXTRABPSIDE"]		# Buffer extra base pairs cut next to the start and end of the BLAST hits
VICINITY = config["VICINITY"]		# Max distance between 5 and 3 end hits to form a haplotype 
IDENTITY = config["IDENTITY"]		# Minimum percentage of identity of BLAST hit to be considered (default 0)

# -------------------------------------------------

# ----------
# Rules not submitted to a job (only needed if you run this in a cluster)
localrules: subset_fasta, subset_gff, gggenomes
# ----------


rule all:
	input:
		"results/SDsynteny.png"

rule subset_fasta:
	input:
		fa = ASSEMBLY
	output:
		fa = "data/haplotypes.fa"
	run:
		with open(output.fa, "w") as out_handle:
			for record in SeqIO.parse(input.fa, "fasta"):
				if record.id in SDscaffolds:
					SeqIO.write(record, out_handle, "fasta")

rule subset_gff:
	""" Make a reduced gff3 for just the SD scaffolds """
	input:
		gff = satDNAgff
	output:
		gff = "data/SDscaffolds.gff3"
	run:
		with open(output.gff, "w") as ofile, open(input.gff, "r") as opengff:
			ofile.write('##gff-version 3\n')

			for line in opengff:
				for scf in SDscaffolds:
					if scf in line:
						ofile.write(line)


# ------- Make alignments --------

rule minimap2:
	input:
		"data/haplotypes.fa"
	output:
		"minimap/haplotypes.paf"
	conda:
		"envs/minimap2.yaml"
	threads: 8,
	# resources: # only needed if you run this in a cluster
	# 	threads = lambda wildcards, threads: threads,
	# 	mem_mb = lambda wildcards, threads: 1500 * threads, 
	# 	time = "2:00:00",
	shell:
		"minimap2 -X -N 50 -p 0.1 -c -B 4 -t {threads} {input} {input} > {output}"
		# "minimap2 -c -B5 -O6 -E3 --rev-only {input} {input} > {output};"
		# -X           skip self and dual mappings (for the all-vs-all mode)
		# -N INT       retain at most INT secondary alignments [5]
		# -p FLOAT     min secondary-to-primary score ratio [0.8]
		# -c           output CIGAR in PAF
		# -B INT       mismatch penalty (larger value for lower divergence) [4]
		# -O INT[,INT] gap open penalty [4,24]
		# -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]

rule reduce_paf:
	"""Modify the paf to have the same number of columns for all entries"""
	input:
		"minimap/haplotypes.paf"
	output:
		"minimap/haplotypes_short.paf"
	shell:
		"cut -f1-12 {input} > {output}"

rule gggenomes: 
	input:
		fa = "data/haplotypes.fa",
		paf = "minimap/haplotypes.paf",
		satDNA = "data/SDscaffolds.gff3",
		genes = genegff
	output:
		plot = "results/SDsynteny.png"
	conda:
		"envs/gggenomes.yaml"	
	script:
		gggenomes

