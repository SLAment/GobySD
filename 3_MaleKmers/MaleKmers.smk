# -*- snakemake -*-

### MaleKmers.smk: Finding male-specific kmers in a Goby genome assembly
############################################################################

# This pipeline uses the kpool output from this command:
# kpool merge -m {table.males} -f {table.females} -o male_specific.tbl

# ==================================================
# S. Lorena Ament-Velasquez
# 2025-03-19
# +++++++++++++++++++++++++++++++++++++++++++++++++

from itertools import combinations

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"

malekmers = config["malekmers"]
SampleIDs = config["SampleIDs"]
Assemblies = config["Assemblies"]
path2genomes = config["path2genomes"]

# Scripts
MaleKmers = config["MaleKmers"]
KmerEnrichment = config["KmerEnrichment"]
# -------------------------------------------------

# ----------
wildcard_constraints:
	sample = '[a-zA-Z0-9_]+'
# ----------

# ----------
# Rules not submitted to a job
localrules: kpool_to_fasta, 
			link_assemblies, 
			indexbam, 
			indexbam_unique, 
			table_hits, 
			collapse_overlapping_kmers, 
			genomebed, 
			makewindows,
			BEDtools_coverage,
			PlotMakeKmers,
			PlotKmerEnrichment
# ----------

# Make dictionary of assemblies
assemblies = dict(zip(SampleIDs, Assemblies))

rule all:
	input:
		"results/FigS7_EnrichmentKmer.png",
		# expand("results/{sample}_scatterplots_n30.png", sample = SampleIDs),
		# expand("results/{sample}_chrpaint_n30.png", sample = SampleIDs),
		# expand("mapping/male_kmers_unique-{sample}.bam.bai", sample = SampleIDs) # to see it in IGV


# -------- Prepare data --------

def get_assembly_name(wildcards):
	return f"{path2genomes}/{assemblies[wildcards.sample]}"

rule link_assemblies:
	""" Give friendlier names to the assemblies """
	input:
		get_assembly_name
	output:
		"data/assemblies/{sample}.fa"
	shell:
		"ln -s {input} {output}"


rule filter_kmers:
	""" Remove some extremes that are probably TEs """
	input:
		kpool = malekmers
	output:
		kpool = "data/male_kmers.tbl"
	params:
		maxcov = 120
	run:
		kmercount = 1
		with open(input.kpool, 'r') as kmers, open(output.kpool, 'w') as newkmers:
			for line in kmers:
				kstring, mcount, fcount = line.rstrip("\n").split("\t")	
				if int(mcount) < params.maxcov and int(fcount) == 0:
					kmerID = f"kmer{kmercount:08d}"
					newline = f"{kmerID}\t{kstring}\t{mcount}\n" # The female counts are irrelevant if they are all 0
					newkmers.write(newline)
				kmercount += 1

rule kpool_to_fasta:
	""" Transform the kmers table from kpool into a fasta file """
	input:
		kpool = "data/male_kmers.tbl"
	output:
		fasta = "data/male_kmers.fa"
	run:
		with open(input.kpool, 'r') as kmers, open(output.fasta, 'w') as fasta:
			for line in kmers: 
				kmerID, kstring, mcount = line.rstrip("\n").split("\t")
				fastaline = f">{kmerID}\n{kstring}\n"
				# fastaline = f">kmer{kmerID:08d}\n{tabs[0]}\n"
				fasta.write(fastaline)

rule indexbwa:
	""" Index genome with BWA """
	input:
		genome = "data/assemblies/{sample}.fa"
	output:
		index = "data/assemblies/{sample}.fa.bwt"
	conda:
		"envs/bwa.yaml"
	threads: 1,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "20GB",
		time = "3:30:00",
	shell:
		"""
		bwa index {input.genome}
		"""

rule bwa_mem:
	input:
		genome = "data/assemblies/{sample}.fa",
		index = "data/assemblies/{sample}.fa.bwt",
		kmers = "data/male_kmers.fa"
	output:
		bam = "mapping/male_kmers_{sample}.bam",
		sai = temp("mapping/male_kmers_{sample}.sai"),
	log:
		"logs/bwa_mem/{sample}.log"
	conda:
		"envs/bwa.yaml"
	threads: 10,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 1500 * threads, #The thin nodes have 128 physical cores and 227 available GB, so just 1.77 GB per core
		time = "2:00:00",
	shell:
		"bwa aln {input.genome} {input.kmers} > mapping/male_kmers_{wildcards.sample}.sai; "
		"(bwa samse {input.genome} mapping/male_kmers_{wildcards.sample}.sai {input.kmers} | samtools view -Su - | samtools sort -O bam -T {wildcards.sample} -@ {threads} > {output.bam}) 2> {log}"
		# "(bwa mem {input.genome} {input.kmers} -t {threads} | samtools view -Su - | samtools sort -O bam -T {wildcards.sample} -@ {threads} > {output.bam}) 2> {log}"

rule indexbam:
	""" Index the BAM file for IGV """
	input:
		"mapping/male_kmers_{sample}.bam"
	output:
		"mapping/male_kmers_{sample}.bam.bai"
	conda:
		"envs/bwa.yaml"
	shell:
		"samtools index {input}"

rule filter_bam:
	input:
		"mapping/male_kmers_{sample}.bam"
	output:
		"mapping/male_kmers_unique-{sample}.bam"
	conda:
		"envs/bwa.yaml"
	shell:
		"samtools view -h {input} | awk '$1 ~ /^@/ || ($5 >= 30 && $0 !~ /NM:i:[1-9]/)' | samtools view -b -o {output}"
		# "samtools view -h {input} | grep -E '@|NM:i:0' | awk '$5 >= 30' | samtools view -b -o {output}"
		# Mismatches are recorded in the NM tag (edit distance), so NM:i:0 means no mismatches.
		# Reads with multiple mappings typically have a MAPQ (mapping quality) value of 0 or low values.

rule indexbam_unique:
	""" Index the BAM file for IGV """
	input:
		"mapping/male_kmers_unique-{sample}.bam",
	output:
		"mapping/male_kmers_unique-{sample}.bam.bai"
	conda:
		"envs/bwa.yaml"
	shell:
		"samtools index {input}"

rule table_hits:
	input:
		"mapping/male_kmers_unique-{sample}.bam",
	output:
		"reports/male_kmers_unique-{sample}.txt"
	conda:
		"envs/bwa.yaml"
	shell:
		"samtools view {input} | awk -v OFS='\\t' '{{print $3, $4, $4+21, $1}}' > {output}"

# ----
def remove_overlap(ranges):
	""" Simplify a list of ranges; I got it from https://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap """
	result = []
	current_start = -1
	current_stop = -1 

	for start, stop in sorted(ranges):
		if start > current_stop:
			# this segment starts after the last segment stops
			# just add a new segment
			result.append( (start, stop) )
			current_start, current_stop = start, stop
		else:
			# current_start already guaranteed to be lower
			current_stop = max(current_stop, stop)
			# segments overlap, replace
			result[-1] = (current_start, current_stop) # SLAV: I modified this to update the stop too.
	return(result)
# ----

rule collapse_overlapping_kmers:
	""" Add flaking base pairs to the bed files to see if that compensates for the missing data """
	input:
		bed = "reports/male_kmers_unique-{sample}.txt"
	output:
		bed = "reports/male_kmers_unique-merged_{sample}.bed"
	run:
		tabs = [line.rstrip("\n").split("\t") for line in open(input.bed)]

		# Make a dictionary of ranges
		beddic = {} # key: chromosome, value: (start, end)

		for tab in tabs:
			contig, start, end, kmer = tab # The kmer info gets lost next
			start = int(start)
			end = int(end)

			if contig in beddic.keys():
				beddic[contig].append([start, end])
			else:
				beddic[contig] = [[start, end]]

		# Reduce the overlaps
		for ctg in beddic.keys():
			beddic[ctg] = remove_overlap(beddic[ctg])

		# Print them in a new file
		with open(output.bed, 'w') as result:
			result.write('#Contig\tStart\tEnd\n') # header
			for ctg in beddic.keys():
				for interval in beddic[ctg]:
					result.write(f'{ctg}\t{interval[0]}\t{interval[1]}\n')

rule genomebed:
	""" Create and index-like bed file for the reference genome """
	input:
		"data/assemblies/{sample}.fa"
	output:
		"tracks/{sample}.bed"
	# conda:
	# 	"envs/bwa.yaml"
	shell:
		"""
		samtools faidx {input}

		cat {input}.fai | awk {{'print $1"\\t"1"\\t"$2'}} > {output} # It has to be tabs or BEDtools won't like it.
		rm {input}.fai 
		"""

rule makewindows:
	""" Use the BEDtools makewindows to get windows of the reference genome """
	input:
		"tracks/{sample}.bed"
	output:
		"tracks/wins_{sample}.bed"
	shell:
		"bedtools makewindows -b {input} -w 50000 -s 50000 > {output}"
		# -w <window_size>
		# -s <step_size>


rule BEDtools_coverage:
	""" Use BEDtools coverage to produce a distribution of male k-mers along the genome """
	# https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
	input:
		chrs = "tracks/wins_{sample}.bed",
		tes = "reports/male_kmers_unique-merged_{sample}.bed"
	output:
		"reports/male_coverage_{sample}.txt"
	shell:
		"""
		bedtools coverage -a {input.chrs} -b {input.tes} | awk '{{ print $1,$2,$3,$7 }}' > {output}
		"""

rule PlotMakeKmers:
	input:
		kmers = "reports/male_coverage_{sample}.txt"
	output:
		points = "results/{sample}_scatterplots_n30.png",
		paint = "results/{sample}_chrpaint_n30.png"
	params:
		minN = 30
	conda:
		"envs/plot.yaml"
	script:
		MaleKmers

rule PlotKmerEnrichment:
	input:
		hap1 = "reports/male_coverage_hap1.txt",
		hap2 = "reports/male_coverage_hap2.txt"
	output:
		plot = "results/FigS7_EnrichmentKmer.png",
	conda:
		"envs/fancyplot.yaml"
	script:
		KmerEnrichment

