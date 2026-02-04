# -*- snakemake -*-

### GobyVariantCalling.smk: Variant calling with Illumina data of gobies to estimate heterozygosity
############################################################################

# ==================================================
# S. Lorena Ament-Velasquez
# 2025-04-22
# +++++++++++++++++++++++++++++++++++++++++++++++++

from Bio import SeqIO

# -------------------------------------------------
# DATA from the configuration file
configfile: "config/config.yaml"
# profilefile: "/cfs/klemming/projects/supr/snic2020-6-175/lore/goby/profile/config.v8+.yaml"

# Sample IDs
SampleIDs = config["SampleIDs"]
# Illumina reads path
path2BAMs = config["path2BAMs"]
# The reference genome
REFGenome = config["REFGenome"]
# The focal scaffold names
targetscfs = config["targetscfs"]
# SD scaffold
SDscfs = config["SDscfs"]
# satDNA gff
satDNA = config["satDNA"]
# TEs gff
TEgff = config["TEgff"]
# Population file
popfile = config["popfile"]
# Pixy window size
WINSIZE = config["WINSIZE"]
# Ploting script
PixyGobySD = config["PixyGobySD"]
# -------------------------------------------------

# ----------
wildcard_constraints:
	scaffold = '[a-zA-Z0-9]+' # to exclude '_'
# ----------


# ----------
# Rules not submitted to a job
localrules: subset_ref, 
			indexsanddict, 
			get_invariants, 
			get_variants, 
			merge_repeats,
			removeoverlap,	
			bgzip_tabix, 
			sample_pop_file, 
			merge_pixy_pi,
			merge_pixy_dxy,
			plot_PixyGobySD

			# sort_gff,
			# pixy_sexes, 
			# get_variants_again, 
			# pixy_pi
# ----------

WINSIZEkb = int(round(WINSIZE/1000))

rule all:
	input:
		"results/Fig5_Heterozygosity_SDscf.png"

# ---------------------------------

rule subset_ref:
	""" Extract the scaffold sequence to use as reference """
	input:
		fa = REFGenome
	output:
		fa = "data/scaffolds/{scaffold}.fa"
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
		"samtools view -b {input} {wildcards.scaffold}.1 > {output}; "
		"samtools index {output}"

# -------------- GATK4 -----------------
rule indexsanddict:
	""" Index reference for GATK """ 
	input:
		genome = "data/scaffolds/{scaffold}.fa",
	output:
		indexsamtools = "data/scaffolds/{scaffold}.fa.fai",
		diction = "data/scaffolds/{scaffold}.dict"
	shell:
		"""
		# Make a reference index
		samtools faidx {input.genome}
		# Make a reference dictionary
		java -jar $PICARD_ROOT/picard.jar CreateSequenceDictionary -R {input.genome} -O {output.diction}
		"""	

rule HaplotypeCaller:
	""" Produce a GVCF file from BAM """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
	input:
		bam = "BAMs/{sample}_{scaffold}.bam",
		ref = "data/scaffolds/{scaffold}.fa",
		indexsamtools = "data/scaffolds/{scaffold}.fa.fai",
		diction = "data/scaffolds/{scaffold}.dict"
	output:
		gvcf = "gvcfs/{sample}_{scaffold}.g.vcf", # I should have used `g.vcf.gz` so it would be compressed by GATK...
	threads: 1,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "60GB", #lambda wildcards, threads: 1200 * threads, #The thin nodes have 128 physical cores and 227 available GB, so just 1.77 GB per core
		time = "6:00:00",
	params:
		JavaMem = "60"
	shell:
		"""
		gatk --java-options "-Xmx{params.JavaMem}G" HaplotypeCaller \\
		-I {input.bam} -R {input.ref} \\
		-O {output.gvcf} \\
		-ERC GVCF -L {wildcards.scaffold}.1

		# --bam-output,-bamout:String   File to which assembled haplotypes should be written  Default value: null.
		# --annotateNDA 	Annotate number of alleles observed
		# --useNewAFCalculator	Use new AF model instead of the so-called exact model
		# --emitRefConfidence GVCF 	Mode for emitting reference confidence scores
		"""

rule makeintervals:
	""" Make an interval file for GenomicsDBImport """
	input:
		"data/scaffolds/{scaffold}.fa.fai"
	output:
		"data/scaffolds/{scaffold}.intervals"
	shell:
		"cat {input} | cut -f1 > {output}"

rule GenomicsDBImport:
	""" Import single-sample GVCFs into GenomicsDB before joint genotyping """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.1/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.php
	# https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport
	input:
		intervals = "data/scaffolds/{scaffold}.intervals",
		gvcfs = expand("gvcfs/{sample}_{{scaffold}}.g.vcf", sample = SampleIDs),
		indexsamtools = "data/scaffolds/{scaffold}.fa.fai",
	output:
		directory("genomicsdb_{scaffold}")
	threads: 1,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "12GB", #lambda wildcards, threads: 1200 * threads, #The thin nodes have 128 physical cores and 227 available GB, so just 1.77 GB per core
		time = "2:00:00",
	params:
		JavaMem = 12
	run:
		# Create a string in the format --variant path/to/gvcf/sample1 --variant path/to/gvcf/sample2 etc...
		variantlist = ""
		for sample in input.gvcfs:
			variantlist += "--variant " + sample + " "
		
		shell("mkdir -p temp")
		# Notice GATK will create the output directory
		gatkcommand = f'gatk GenomicsDBImport --java-options "-Xmx{params.JavaMem}G" --tmp-dir temp --genomicsdb-workspace-path {output[0]} -L {input.intervals} {variantlist}'
		# gatkcommand = f'gatk GenomicsDBImport --java-options "-Xmx{params.JavaMem}G" --tmp-dir temp --genomicsdb-workspace-path {output[0]} -L {input.intervals} {variantlist} --sequence-dictionary {input.diction}'
		shell(gatkcommand) # execute

rule GenotypeGVCFs_allsites:
	""" Perform joint genotyping and produce an all-sites vcf file """
	# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php
	# https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport
	# https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html
	input:
		my_database = "genomicsdb_{scaffold}",
		ref = "data/scaffolds/{scaffold}.fa",
		diction = "data/scaffolds/{scaffold}.dict"
	output:
		rawvcf = "vcfs/Raw_{scaffold}.vcf"
	threads: 1, # it doesn't matter, it's just the memory that matters
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = "40GB", #lambda wildcards, threads: 1200 * threads, #The thin nodes have 128 physical cores and 227 available GB, so just 1.77 GB per core
		time = "2:00:00",
	params:
		JavaMem = 40
	shell:
		"""
		gatk --java-options "-Xmx{params.JavaMem}G" GenotypeGVCFs \\
		-R {input.ref} \\
		-V gendb://{input.my_database} \\
		-O {output.rawvcf} \\
		-all-sites
		"""
# --create-output-variant-index	If true, create a VCF index when writing a coordinate-sorted VCF file.
# --use-new-qual-calculator / -new-qual 	Use the new AF model instead of the so-called exact model. Default: true
# By default, GATK HaplotypeCaller and GenotypeGVCFs do not emit variants with QUAL < 10, controlled with -stand-call-conf

# Nima added 
# -G StandardAnnotation
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622
# https://gatk.broadinstitute.org/hc/en-us/articles/360050814612-HaplotypeCaller
# --annotation-group, -G:String  One or more groups of annotations to apply to variant calls. This argument may be specified 0 or more times. Default value: null. Possible Values: {AlleleSpecificAnnotation, AS_StandardAnnotation, ReducibleAnnotation, StandardAnnotation, StandardHCAnnotation, StandardMutectAnnotation}

rule fix_missingdata: # I used 4.5.0.0, but it's fixed in 4.6.0.0
	""" GATK changed the formatting of missing data and now it's not standard, ugh. Fix it. """
	# https://github.com/broadinstitute/gatk/issues/8328
	# https://twitter.com/ksamuk/status/1754698564884996590
	# https://samtools.github.io/bcftools/bcftools.html
	input:
		"vcfs/Raw_{scaffold}.vcf"
	output:
		"vcfs/Raw_{scaffold}_fixed.vcf.gz"
	shell:
		# "bcftools +setGT {input} -- -t q -n . -i 'FORMAT/DP=0 | SMPL_MAX(FORMAT/PL)=0' > {output}"
		"bcftools +setGT {input} -- -t q -n . -e 'FMT/DP>=1' | bgzip -c > {output}"
	# -t, --target-gt <type>      Genotypes to change
	#	q    .. select genotypes using -i/-e options
	# -n, --new-gt <type>         Genotypes to set
	# 	.    .. partially or completely missing
	# -e, --exclude <expr>        Exclude a genotype if true (requires -t q)
	# when prefixed with SMPL_ (or "s" for brevity, e.g. SMPL_MAX or sMAX), they will evaluate to a vector of per-sample values when applied on FORMAT tags:
	# -i, --include EXPRESSION include only sites for which EXPRESSION is true.

rule get_invariants:
	input:
		"vcfs/Raw_{scaffold}_fixed.vcf.gz"
	output:
		"vcfs/Raw_{scaffold}_fixed_invar.vcf.gz"
	shell:
		"vcftools --gzvcf {input} --max-maf 0 --min-meanDP 20 --max-meanDP 300 --recode --stdout | bgzip -c > {output}; "
		"tabix {output}" # necessary for merging

rule get_variants:
	input:
		"vcfs/Raw_{scaffold}_fixed.vcf.gz"
	output:
		"vcfs/Raw_{scaffold}_fixed_var.vcf.gz"
	shell:
		"vcftools --gzvcf {input} --mac 1 --min-meanDP 20 --max-meanDP 300 --recode --stdout | bgzip -c > {output}; "
		"bcftools index -t {output}" # it will be necessary for GATK (it could have been with tabix)
		# --mac	Minor Allele Count

rule MarkFilter:
	""" Mark variant sites that pass the filter """
	# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
	input:
		vcf = "vcfs/Raw_{scaffold}_fixed_var.vcf.gz"
	output:
		vcf = "vcfs/Raw_{scaffold}_fixed_var_marked.vcf"
	shell:
		"""
		gatk VariantFiltration \\
		--variant {input.vcf} \\
		-O {output.vcf} \\
		-filter "QD < 2.0" --filter-name "QD2" \\
		-filter "FS > 10.0" --filter-name "FS10" \\
		-filter "ReadPosRankSum < -3.0" --filter-name "ReadPosRankSum-3" \\
		-filter "MQRankSum  < -6.0 " --filter-name "MQRankSum-6" \\
		-filter "SOR > 3.0" --filter-name "SOR-3" \\
		-filter "MQ < 40.0" --filter-name "MQ-40" 
		"""

rule getPASSsites:
	""" Remove sites rejected in VariantFiltration and put them together """
	input:
		"vcfs/Raw_{scaffold}_fixed_var_marked.vcf",
	output:
		"vcfs/Raw_{scaffold}_fixed_var_pass.vcf",
	shell:
		"""
		bcftools view -h {input} > {output} # header
		bcftools view -H {input} | awk '($7 == "PASS")' >> {output} # header
		"""

rule compress_vars:
	input:
		vcf = "vcfs/Raw_{scaffold}_fixed_var_pass.vcf",
	output:
		vcf = "vcfs/Raw_{scaffold}_fixed_var_pass.vcf.gz"
	shell:
		"bgzip {input.vcf} && tabix -p vcf {output.vcf}"


rule merge_vcfs:
	""" Combine the two VCFs using bcftools concat """
	input:
		"vcfs/Raw_{scaffold}_fixed_var_pass.vcf.gz",
		"vcfs/Raw_{scaffold}_fixed_invar.vcf.gz"
	output:
		"vcfs/{scaffold}_qc.vcf.gz"
	shell:
		"bcftools concat --allow-overlaps {input} -O z -o {output}"

# -------------- Remove satDNA -----------------

rule merge_repeats:
	input:
		sat = satDNA,
		TEs = TEgff
	output:
		temp("temp/repeats.gff")
	shell:
		"cat {input} | grep -v '^##' >> {output}"

# rule sort_gff: # not necessary in the end
# 	input:
# 		"temp/repeats.gff"
# 	output:
# 		temp("temp/repeats_sorted.gff")
# 	shell:
# 		"sort -k1,1 -k2,2n {input} > {output}"

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

rule removeoverlap:
	""" Fuse overlapping ranges in a BED file """
	# Basically the same as the script totalcovergff.py but for bed files
	input:
		gff = "temp/repeats.gff"
	output:
		bed = "reports/repeats_merged.bed"
	run:
		# Make a dictionary of ranges
		beddic = {} # key: chromosome, value: (start, end)

		with open(input.gff, 'r') as file: 
			for line in file:
				if '##' in line: # although not necessary here
					pass
				elif line not in ['\n', '\r\n']: # Ignore empty lines
					cols = line.rstrip("\n").split("\t")

					contig = cols[0]
					start = int(cols[3])
					end = int(cols[4])

					if contig in list(beddic.keys()):
						beddic[contig].append([start, end])
					else: # contig is new
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


rule bedtools_repeats:
	""" Filter vcf with the repeats from RepeatMasker """
	input:
		vcf = "vcfs/{scaffold}_qc.vcf.gz",
		gfffile = "reports/repeats_merged.bed"
	output:
		filteredvcf = "vcfs/{scaffold}_qc-nosatDNA.vcf",
	shell:
		"""
		# Get the header
		bcftools view -h {input.vcf} > {output.filteredvcf}
		
		# Filter out the repeats
		bedtools intersect -a {input.vcf} -b {input.gfffile} -v >> {output.filteredvcf}
		"""

rule bgzip_tabix:
	""" Compress and index vcf file """
	input:
		"vcfs/{scaffold}_qc-nosatDNA.vcf"
	output:
		"vcfs/{scaffold}_qc-nosatDNA.vcf.gz"
	shell:
		"bgzip {input}; "
		"tabix -p vcf {output}" # Necessary for pixy

# -------------- pixy -----------------

rule pixy_sexes: # it's super fast!
	""" Calculate divergence between females and males """
	input:
		vcf = "vcfs/{scaffold}_qc-nosatDNA.vcf.gz",
		popu = popfile
	output:
		"pixy/{scaffold}_sex_{WINSIZEkb}kb_dxy.txt"
	conda:
		"pixyenv" # https://stackoverflow.com/questions/59107413/activating-existing-conda-enviornments-in-snakemake
		# "envs/pixy.yaml"
	params:
		winsize = WINSIZE
	threads: 8,
	# resources:
	# 	threads = lambda wildcards, threads: threads,
	# 	mem_mb = lambda wildcards, threads: 1200 * threads, #The thin nodes have 128 physical cores and 227 available GB, so just 1.77 GB per core
	# 	time = "1:00:00",
	shell:
		"""
		pixy --stats dxy \\
		--vcf {input.vcf} \\
		--populations {input.popu} \\
		--window_size {params.winsize} \\
		--n_cores {threads} \\
		--output_folder pixy \\
		--output_prefix {wildcards.scaffold}_sex_{wildcards.WINSIZEkb}kb \\
		--chromosomes {wildcards.scaffold}.1
		"""
		# pixy ignores non-biallelic sites and INDELs,

rule sample_pop_file:
	""" Prepare a population file to calculate heterozygosity per sample """
	output:
		popu = "temp/samples_popfile.txt"
	run:
		with open(output.popu, "w") as ofile:
			for sample in SampleIDs:
				newline = f"{sample}\t{sample}\n"
				ofile.write(newline)

# rule get_variants_again:
# 	""" Get the variant sites again after the filtering """
# 	input:
# 		"vcfs/{scaffold}_qc-nosatDNA.vcf.gz"
# 	output:
# 		"vcfs/{scaffold}_qc-nosatDNA_var.vcf.gz"
# 	shell:
# 		"vcftools --gzvcf {input} --mac 1 --recode --stdout | bgzip -c > {output}; "
# 		"bcftools index -t {output}" # it will be necessary for GATK (it could have been with tabix)
# 		# --mac	Minor Allele Count

rule pixy_pi: # Takes slightly longer
	""" Calcualte heterozygosity of each sample """
	input:
		vcf = "vcfs/{scaffold}_qc-nosatDNA.vcf.gz",
		# vcf = "vcfs/{scaffold}_qc-nosatDNA_var.vcf.gz",
		popu = "temp/samples_popfile.txt"
	output:
		"pixy/{scaffold}_{WINSIZEkb}kb_pi.txt"
	conda:
		"pixyenv" # https://stackoverflow.com/questions/59107413/activating-existing-conda-enviornments-in-snakemake
		# "envs/pixy.yaml"
	params:
		winsize = WINSIZE
	threads: 8,
	resources:
		threads = lambda wildcards, threads: threads,
		mem_mb = lambda wildcards, threads: 1200 * threads, #The thin nodes have 128 physical cores and 227 available GB, so just 1.77 GB per core
		time = "1:00:00",
	shell:
		"""
		pixy --stats pi \\
		--vcf {input.vcf} \\
		--populations {input.popu} \\
		--window_size {params.winsize} \\
		--n_cores {threads} \\
		--output_folder pixy \\
		--output_prefix {wildcards.scaffold}_{wildcards.WINSIZEkb}kb \\
		--chromosomes {wildcards.scaffold}.1 #--bypass_invariant_check
		"""
		# pixy ignores non-biallelic sites and INDELs

rule merge_pixy_pi:
	input:
		expand("pixy/{scaffold}_{{WINSIZEkb}}kb_pi.txt", scaffold = targetscfs),
	output:
		"reports/{WINSIZEkb}kb_pi.txt"
	shell:
		"printf 'pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing\n' > {output}; "
		"cat {input} | grep -v 'indow_pos_1' >> {output}"

rule merge_pixy_dxy:
	input:
		expand("pixy/{scaffold}_sex_{{WINSIZEkb}}kb_dxy.txt", scaffold = targetscfs)
	output:
		"reports/sex_{WINSIZEkb}kb_dxy.txt"
	shell:
		"printf 'pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_dx\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing\n' > {output}; "
		"cat {input} | grep -v 'indow_pos_1' >> {output}"

rule plot_PixyGobySD:
	""" Plot the coverage along the scaffold """
	input:
		TEs = TEgff,
		satDNA = satDNA,
		het = expand("reports/{WINSIZEkb}kb_pi.txt", WINSIZEkb = WINSIZEkb),
		dxy = expand("reports/sex_{WINSIZEkb}kb_dxy.txt", WINSIZEkb = WINSIZEkb),
	output:
		het = "results/Fig5_Heterozygosity_SDscf.png",
		dxy = "results/Divergence_SDscf.png",
	conda:
		"envs/plot.yaml"
	script:
		PixyGobySD	



