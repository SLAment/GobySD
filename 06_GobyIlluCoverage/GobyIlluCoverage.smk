# -*- snakemake -*-

### GobyIlluCoverage: Calculating coverage along the goby's genome
############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# 2025-10-08
# +++++++++++++++++++++++++++++++++++++++++++++++++

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"
profilefile: "../4_GobyMapping/profile/config.v8+.yaml"

# Sample IDs
SampleIDs = config["SampleIDs"]
# Illumina reads path
path2bams = config["path2bams"]
# The reference genome
REFGenome = config["REFGenome"]
# satDNA annotation 
satDNA = config["satDNA"]
# Scripts
GobyIlluCovDist = config["GobyIlluCovDist"]

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: linkref,  
			genomebed, 
			makewindows, 
			makeindex, 
			aggregate_coverage, 
			subset_windows, 
			merge_tables, 
			plot_GobyIlluCovDist
# ----------

rule all:
	input:
		"results/Coverage_density_chrs.png",
		"results/Coverage_scan_chrs.png",
		# expand("scf_bams/{sample}", sample = SampleIDs),

# ---------------------------------

rule linkref:
	""" Make a link to the reference so it doesn't put new stuff in the ref folder """
	input:
		REFGenome
	output:
		"data/genome.fa"
	shell:
		"ln -s {input} {output}"

# -------------- Coverage -----------------

rule genomebed:
	""" Create and index-like bed file for the reference genome """
	input:
		"data/genome.fa"
	output:
		"tracks/genome.bed"
	shell:
		"""
		samtools faidx {input}

		cat {input}.fai | awk {{'print $1"\\t"1"\\t"$2'}} > {output} # It has to be tabs or BEDtools won't like it.
		rm {input}.fai 
		"""

rule makewindows:
	""" Use the BEDtools makewindows to get windows of the reference genome """
	input:
		"tracks/genome.bed"
	output:
		"tracks/wins_genome.bed"
	shell:
		"bedtools makewindows -b {input} -w 30000 -s 30000 > {output}"
		# -w <window_size>
		# -s <step_size>


rule filter_bam:
	""" Remove questionable mappings from the BAM file """
	input:
		bam = path2bams + "/{sample}/{sample}.sorted.debup.bam",
	output:
		bam = "mapping/{sample}/{sample}.sorted.debup.filtered.bam"
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

checkpoint scf2BAM:
	""" Extract each contig as a BAM file because doing it all at a time runs out of memory """
	input:
		bed = "tracks/wins_genome.bed",
		bam = "mapping/{sample}/{sample}.sorted.debup.filtered.bam",
		index = "mapping/{sample}/{sample}.sorted.debup.filtered.bam.bai"
	output:
		directory("scf_bams/{sample}")
	threads: 1,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 1500 * threads,
		time = "5:00:00",
	run:
		tabs = list(set([line.rstrip("\n").split("\t")[0] for line in open(input.bed, 'r') if 'OZ251433.1' not in line])) # ignore the mitochondria

		# Make the expected output of the checkpoint:
		shell(f"mkdir -p scf_bams/{wildcards.sample}")

		for scf in tabs:
			subbamname = f"scf_bams/{wildcards.sample}/{scf}/{scf}.bam"
			shell(f"mkdir -p scf_bams/{wildcards.sample}/{scf}")

			# Subset the BAM file
			shell(f"samtools view -b {input.bam} {scf} > {subbamname}")
			# -b tells it to output in BAM format		
			
			# Make an index too
			shell(f"samtools index {subbamname}")

			# # This idea didn't work because Snakemake remakes the folder as it's the output that it's tracking...
			# if os.path.exists(os.path.join(os.getcwd(), subbamname)): # leave it alone if it already exists
			# 	shell(f"touch {subbamname}") # just fool snakemake to think it's new
			# 	print(f"{subbamname} already exists, leaving it as is ...")

			# else: # it doesn't exist, Subset the BAM file
			# 	shell(f"mkdir -p scf_bams/{wildcards.sample}/{scf}")
			# 	shell(f"samtools view -b {input.bam} {scf} > {subbamname}")
			# 	# -b tells it to output in BAM format			


# https://edwards.flinders.edu.au/how-to-use-snakemake-checkpoints/
# https://evodify.com/snakemake-checkpoint-tutorial/
def get_scf_bams(wildcards):
	checkpoint_output = checkpoints.scf2BAM.get(**wildcards).output[0]  # The output folder's name
	SCF1, SCF2 = glob_wildcards(os.path.join(checkpoint_output, "{scf}/{scf2}.bam")) # SCF1 = SCF2 but the dictionary fails if I treat them the same
	return expand("scf_bams/{sample}/{scf}/{scf}_covwin30kb.bed", scf= SCF1, sample = wildcards.sample)


rule subset_windows:
	input:
		bed = "tracks/wins_genome.bed",
	output:
		"tracks/wins_{scf}.bed"
	shell:
		"grep {wildcards.scf} {input} > {output}"

rule covperwin:
	""" Use BEDtools to calculate average depth of coverage """
	input:
		bed = "tracks/wins_{scf}.bed",
		bam = "scf_bams/{sample}/{scf}/{scf}.bam",
	output:
		"scf_bams/{sample}/{scf}/{scf}_covwin30kb.bed"
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
		get_scf_bams
	output:
		"reports/{sample}_covwin30kb.bed"
	shell:
		"cat {input} > {output}"

rule merge_tables:
	input:
		expand("reports/{sample}_covwin30kb.bed", sample = SampleIDs)
	output:
		report = "results/Goby_covwin30kb.bed"
	run:
		with open(output.report, 'w') as ofile:
			# Header
			ofile.write('seqid\tstart\tend\tCoverage\tSample\tSex\n')
			for file in input:
				with open(file, 'r') as covwin:
					# Infer the sample 
					sample = file.split("/")[1].split('_')[0]
					# Infer the sex
					if sample[-1] == "f":
						sex = "Female"
					elif sample[-1] == "m":
						sex = "Male"

					for line in covwin:
						newline = line.rstrip('\n') + f"\t{sample}\t{sex}\n"
						ofile.write(newline)

rule plot_GobyIlluCovDist:
	""" Plot the coverage along chromosomes """
	input:
		covwin = "results/Goby_covwin30kb.bed",
		gff = satDNA
	output:
		density = "results/Coverage_density_chrs.png",
		chrs = "results/Coverage_scan_chrs.png",
	conda:
		"envs/plot.yaml"
	script:
		GobyIlluCovDist	


