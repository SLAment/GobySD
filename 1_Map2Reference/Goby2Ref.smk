# -*- snakemake -*-

### Goby2Ref
#############################################################################
# Comparing our assembly to the published chr-level genome assembly from ATLASea (fGobFla1, BioProject PRJEB88435)

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/05/19
# ++++++++++++++++++++++++++++++++++++++++++++++

import re
from Bio import SeqIO

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"

# Data
REF = config["reference"]
NAMEREF = config["nameref"]
QUERY = config["query"]
NAMEQUERY = config["namequery"]

# Scripts
dotplot = config["dotplot"]

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: get_DeltaSenseCorrector, DeltaSenseCorrector, get_chrlens, dotplot
# ----------


rule all:
	input:
		f"results/{NAMEQUERY}-vs-{NAMEREF}_Dotplot.png", # Fig. S2
		f"results/{NAMEQUERY}-vs-{NAMEREF}_Dotplot_chr16.png", # Fig. S5
		f"results/{NAMEQUERY}-vs-{NAMEREF}_chr16_bar.png" # part of Fig. 3
		# f"results/{NAMEQUERY}.fa", # The final TH1 assembly

# -------- Align --------

rule mummer:
	""" Align the draft assembly to the chromosome-level assembly of the ref goby """
	input:
		reference = REF,
		query = QUERY
	output:
		delta = "mummer/{NAMEQUERY}-vs-{NAMEREF}.delta",
		deltafilter = "mummer/{NAMEQUERY}-vs-{NAMEREF}.filter",
		coords = "mummer/{NAMEQUERY}-vs-{NAMEREF}.coords",
		coordsfilter = "mummer/{NAMEQUERY}-vs-{NAMEREF}.filter.coords",	
	threads: 3,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "30GB",
		time = "2-00:00:00",
		partition = "long" # job > 24 hours
	shell:
		"""
		echo "MUMmer alignment ..."
		nucmer -b 200 -c 65 -p mummer/{NAMEQUERY}-vs-{NAMEREF} {input.reference} {input.query}

		echo "Filter the delta"
		echo delta-filter -q {output.delta}
		delta-filter -q {output.delta} > {output.deltafilter}

		# To view a summary of all the alignments produced by NUCmer
		echo "Running show-coords"
		# For Ribbon http://genomeribbon.com/
		echo "...for Ribbon"
		show-coords -r -lTH {output.delta} > {output.coords}
		show-coords -r -lTH {output.deltafilter} > {output.coordsfilter}

		"""
		# --mum  Use anchor matches that are unique in both the reference and query
		# --mumreference  Use anchor matches that are unique in in the reference
        #           but not necessarily unique in the query (default behavior)
        # -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
        # -b|breaklen     Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)


rule get_DeltaSenseCorrector:
	""" Get script from my GitHub """
	output:
		"scripts/DeltaSenseCorrector.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/refs/heads/master/Miscellaneous/DeltaSenseCorrector.py"

rule DeltaSenseCorrector:
	input:
		script = "scripts/DeltaSenseCorrector.py",
		delta = f"mummer/{NAMEQUERY}-vs-{NAMEREF}.delta",
		fasta = QUERY
	output:
		fa = "results/{NAMEQUERY}.fa"
	shell:
		"python {input.script} {input.delta} {input.fasta} --verbose -o {output.fa}"

rule mummer_again:
	""" Align the draft assembly to the chromosome-level assembly of the ref goby """
	input:
		reference = REF,
		query = "results/{NAMEQUERY}.fa"
	output:
		delta = "mummer2/{NAMEQUERY}-vs-{NAMEREF}_sc.delta",
		deltafilter = "mummer2/{NAMEQUERY}-vs-{NAMEREF}_sc.filter",
		coords = "mummer2/{NAMEQUERY}-vs-{NAMEREF}_sc.coords",
		coordsfilter = "mummer2/{NAMEQUERY}-vs-{NAMEREF}_sc.filter.coords",	
	threads: 3,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "30GB",
		time = "2-00:00:00",
		partition = "long" # job > 24 hours
	shell:
		"""
		echo "MUMmer alignment ..."
		nucmer -b 200 -c 65 -p mummer2/{NAMEQUERY}-vs-{NAMEREF}_sc {input.reference} {input.query}

		echo "Filter the delta"
		echo delta-filter -q {output.delta}
		delta-filter -q {output.delta} > {output.deltafilter}

		# To view a summary of all the alignments produced by NUCmer
		echo "Running show-coords"
		# For Ribbon http://genomeribbon.com/
		echo "...for Ribbon"
		show-coords -r -lTH {output.delta} > {output.coords}
		show-coords -r -lTH {output.deltafilter} > {output.coordsfilter}

		"""
		# --mum  Use anchor matches that are unique in both the reference and query
		# --mumreference  Use anchor matches that are unique in in the reference
        #           but not necessarily unique in the query (default behavior)
        # -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
        # -b|breaklen     Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)

rule get_chrlens:
	input:
		fasta = REF,
	output:
		chrlens = temp("temp/chrlens.txt")
	run:
		with open(output.chrlens, 'w') as ofile:
			for seq_record in SeqIO.parse(input.fasta, "fasta"):
				ofile.write(f'{seq_record.id}\t{len(seq_record)}\n')

rule dotplot:
	""" Make dotplot of the REF vs QUERY """
	input:
		coords = "mummer2/{NAMEQUERY}-vs-{NAMEREF}_sc.filter.coords",
		chrlens = "temp/chrlens.txt",
		script = dotplot
	output:
		big = "results/{NAMEQUERY}-vs-{NAMEREF}_Dotplot.png",
		chr16 = "results/{NAMEQUERY}-vs-{NAMEREF}_Dotplot_chr16.png",
		chr16bar = "results/{NAMEQUERY}-vs-{NAMEREF}_chr16_bar.png" # part of Fig. 3
	conda:
		"envs/dotplots.yaml"	
	script:
		dotplot	
