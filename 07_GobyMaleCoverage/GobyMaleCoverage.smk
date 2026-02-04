# -*- snakemake -*-

### GobyMaleCoverage: Calculating coverage along the goby's genome
############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# 2025-04-09
# +++++++++++++++++++++++++++++++++++++++++++++++++

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"
profilefile: "../4_GobyMapping/profile/config.v8+.yaml"

# Sample IDs
SampleIDs = config["SampleIDs"]
# Illumina reads path
path2BAMs = config["path2BAMs"]
# The reference genome
REFGenome = config["REFGenome"]
# List of scaffolds to calculate coverage (because using them all is too much)
focalscf = config["focalscf"]
# Gene Annotation
genegff = config["genegff"]
# satDNA annotation
satDNA = config["satDNA"]

# Window size
WINSIZE = config["WINSIZE"]
# Scripts
GobyIlluCov = config["GobyIlluCov"]
SDstart = config["SDstart"]

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: subset_ref, genomebed, makewindows, makeindex, aggregate_coverage, Plot_GobyIlluCov
# ----------

# ----------
wildcard_constraints:
	scaffold = '[a-zA-Z0-9.]+' # to exclude '_'
# ----------

from Bio import SeqIO

WINSIZEkb = int(round(WINSIZE/1000))

rule all:
	input:
		expand("results/{scaffold}_covwin{WINSIZEkb}kb_coverage_scan.png", scaffold = focalscf, WINSIZEkb = WINSIZEkb),

# ---------------------------------

rule subset_ref:
	""" Extract the scaffold sequence to use as reference """
	input:
		fa = REFGenome
	output:
		fa = temp("data/scaffolds/{scaffold}.fa")
	run:
		with open(output.fa, "w") as ofile:
			for record in SeqIO.parse(input.fa, "fasta"):
				if wildcards.scaffold in record.id:
					SeqIO.write(record, ofile, "fasta")

rule subset_BAM:
	input:
		path2BAMs + "/{sample}/{sample}.sorted.debup.filtered.bam"
	output:
		"BAMs/{sample}_{scaffold}.bam"
	shell:
		"samtools view -b {input} {wildcards.scaffold} > {output}; "
		"samtools index {output}"


# -------------- Coverage -----------------

rule genomebed:
	""" Create and index-like bed file for the reference genome """
	input:
		"data/scaffolds/{scaffold}.fa"
	output:
		"tracks/{scaffold}.bed"
	shell:
		"""
		samtools faidx {input}

		cat {input}.fai | awk {{'print $1"\\t"1"\\t"$2'}} > {output} # It has to be tabs or BEDtools won't like it.
		rm {input}.fai 
		"""

rule makewindows:
	""" Use the BEDtools makewindows to get windows of the reference genome """
	input:
		"tracks/{scaffold}.bed"
	output:
		"tracks/{scaffold}_wins{WINSIZEkb}kb.bed"
	shell:
		"bedtools makewindows -b {input} -w {WINSIZE} -s {WINSIZE} > {output}"
		# -w <window_size>
		# -s <step_size>


rule filter_bam:
	""" Remove questionable mappings from the BAM file """
	input:
		bam = path2BAMs + "/{sample}/{sample}.sorted.debup.bam",
	output:
		bam = temp("mapping/{sample}/{sample}.sorted.debup.filtered.bam")
	threads: 2,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 1500 * threads,
		time = "3:00:00",
	shell:
		"samtools view -b -F 0x904 -F 0x400 {input.bam} > {output.bam}"
		# -F 0x904: removes secondary (0x100), supplementary (0x800), and unmapped (0x4) reads. (supplementary: different parts of a read to different locations)
		# -F 0x400: removes duplicates.
		# -q 20: removes reads with mapping quality < 20.

rule makeindex:
	input:
		"mapping/{sample}/{sample}.sorted.debup.filtered.bam"
	output:
		"mapping/{sample}/{sample}.sorted.debup.filtered.bam.bai"
	shell:
		"samtools index {input}"

rule scf2BAM:
	input:
		bam = "mapping/{sample}/{sample}.sorted.debup.filtered.bam",
		index = "mapping/{sample}/{sample}.sorted.debup.filtered.bam.bai"		
	output:
		"scf_bams/{sample}/{scaffold}_filtered.bam"
	threads: 1,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 1500 * threads,
		time = "2:00:00",
	shell:
		"samtools view -b {input.bam} {wildcards.scaffold} > {output}; "
		"samtools index {output}"

rule covperwin:
	""" Use BEDtools to calculate average depth of coverage """
	input:
		bed = "tracks/{scaffold}_wins{WINSIZEkb}kb.bed",
		bam = "scf_bams/{sample}/{scaffold}_filtered.bam"
	output:
		"BEDs/{scaffold}_{sample}_covwin{WINSIZEkb}kb.bed"
	threads: 2,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "60GB",
		time = "1:30:00",
	shell:
		"bedtools coverage -a {input.bed} -b {input.bam} -mean > {output}"

rule aggregate_coverage:
	""" Bring together the coverage estimates """
	input:
		expand("BEDs/{{scaffold}}_{sample}_covwin{{WINSIZEkb}}kb.bed", sample = SampleIDs)
	output:
		report = "reports/{scaffold}_covwin{WINSIZEkb}kb.bed"
	run:
		with open(output.report, 'w') as ofile:
			# Header
			ofile.write('seqid\tstart\tend\tCoverage\tSample\tSex\n')
			for file in input:
				with open(file, 'r') as covwin:
					# Infer the sample 
					sample = file.split("/")[1].split('_')[1]
					# Infer the sex
					if sample[-1] == "f":
						sex = "Female"
					elif sample[-1] == "m":
						sex = "Male"

					for line in covwin:
						newline = line.rstrip('\n') + f"\t{sample}\t{sex}\n"
						ofile.write(newline)

rule plot_GobyIlluCov:
	""" Plot the coverage along the scaffold """
	input:
		covwin = "reports/{scaffold}_covwin{WINSIZEkb}kb.bed",
		gff = satDNA,
		gene = genegff
	output:
		chrs = "results/{scaffold}_covwin{WINSIZEkb}kb_coverage_scan.png",
	conda:
		"envs/plot.yaml"
	params:
		start = SDstart
	script:
		GobyIlluCov	


