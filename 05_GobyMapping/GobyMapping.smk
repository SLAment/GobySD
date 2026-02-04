# -*- snakemake -*-

### GobyMapping.smk: Preparing BAM files for variant calling
############################################################################

# ==================================================
# S. Lorena Ament-Velasquez
# 2025-03-21
# +++++++++++++++++++++++++++++++++++++++++++++++++


# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"

# Sample IDs
SampleIDs = config["SampleIDs"]
# Illumina reads path
path2Illumina = config["path2Illumina"]
# The reference genome
REFGenome = config["REFGenome"]
# All the Goby Illumina was sequenced in the same lane
flowcell = config["flowcell"]

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: linkref
# ----------


rule all:
	input:
		expand("mapping/{sample}/{sample}.sorted.debup.bam", sample = SampleIDs[0]),

# ---------------------------------

rule linkref:
	""" Make a link to the reference so it doesn't put new stuff in the ref folder """
	input:
		REFGenome
	output:
		"data/genome.fa"
	shell:
		"ln -s {input} {output}"

rule indexbwa:
	""" Index genome with BWA """
	input:
		genome = "data/genome.fa"
	output:
		index = "data/genome.fa.bwt"
	threads: 3,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "6GB", #lambda wildcards, threads: 1200 * threads, #The thin nodes have 128 physical cores and 227 available GB, so just 1.77 GB per core
		time = "6:00:00",
	shell:
		"""
		bwa index {input.genome}
		"""

rule bwa_mem:
	""" Map Illumina reads with BWA """
	input:
		genome = "data/genome.fa",
		index = "data/genome.fa.bwt",
		read1 = path2Illumina + "/{sample}/{sample}_postQC.1.fq.gz",
		read2 = path2Illumina + "/{sample}/{sample}_postQC.2.fq.gz",
	output:
		bam = temp("mapping/{sample}/{sample}.bam.sorted"),
	log:
		"logs/bwa_mem/{sample}.log"
	threads: 20,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "80GB", #lambda wildcards, threads: 1200 * threads, #The thin nodes have 128 physical cores and 227 available GB, so just 1.77 GB per core
		time = "3-00:00:00",
	params:
		rg = lambda wildcards: f"@RG\\tID:{flowcell}\\tSM:{wildcards.sample}\\tPL:ILLUMINA", # https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
	shell:
		"""
		(bwa mem {input.genome} {input.read1} {input.read2} -t {threads} -R '{params.rg}' -M | samtools view -Su - | samtools sort -l 5 -O bam -T {wildcards.sample} -@ {threads} > {output.bam}) 2> {log}
		# -l 5 following Doug
		"""

rule markduplicates:
	""" Mark duplicates in BAM """
	# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
	input:
		bam = "mapping/{sample}/{sample}.bam.sorted"
	output:
		mdoutput = "mapping/{sample}/{sample}.sorted.debup.bam",
		mdmetrics = "mapping/{sample}/{sample}.sorted.metrics.txt"
	threads: 1,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "40GB", #lambda wildcards, threads: 1200 * threads, #The thin nodes have 128 physical cores and 227 available GB, so just 1.77 GB per core
		time = "1-00:00:00",
	shell:
		"""
		# Using normal Picard
		java -jar $PICARD_ROOT/picard.jar MarkDuplicates -I {input.bam} -O {output.mdoutput} -M {output.mdmetrics} --ASSUME_SORT_ORDER coordinate --CREATE_INDEX true --REMOVE_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --TMP_DIR $PDC_TMP
		"""	
		# # VALIDATION_STRINGENCY=ValidationStringency
		# #                               Validation stringency for all SAM files read by this program.  Setting stringency to
		# #                               SILENT can improve performance when processing a BAM file in which variable-length data
		# #                               (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This
		# #                               option can be set to 'null' to clear the default value. Possible values: STRICT,
		# #                               LENIENT, SILENT
		# # CREATE_INDEX=Boolean          Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value:
		# #                               false. This option can be set to 'null' to clear the default value. Possible values:
		# #                               true, false
		# # TMP_DIR (File)  Default value: null. This option may be specified 0 or more times.
		# # --REMOVE_DUPLICATES:Boolean   If true do not write duplicates to the output file instead of writing them with
        # # 		                      appropriate flags set.  Default value: false. Possible values: {true, false}
        # # --OPTICAL_DUPLICATE_PIXEL_DISTANCE:Integer
        # #                               The maximum offset between two duplicate clusters in order to consider them optical
        # #                               duplicates. The default is appropriate for unpatterned versions of the Illumina platform.
        # #                               For the patterned flowcell models, 2500 is moreappropriate. For other platforms and
        # #                               models, users should experiment to find what works best.  Default value: 100.
        # #	--ASSUME_SORT_ORDER,-ASO:SortOrder
        # #                               If not null, assume that the input file has this order even if the header says otherwise.
        # #                               Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate,
        # #                               unknown}  Cannot be used in conjunction with argument(s) ASSUME_SORTED (AS)

