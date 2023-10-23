configfile:"sample_list.json"
PATH="Path directory"
REF="ucsc.hg19.fasta"
BED="human_hg19.bed"
KNOWN_VARIANTS="dbsnp_138.hg19.vcf"
KNOWN_INDELS="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
MUTECT_GERMLINE="af-only-gnomad.raw.sites.hg19.vcf"
PE_REF="TruSeq2-PE.fa"
## WES

rule FASTQC:
        input:
                r1=SAMPLE_DIR_ORIGIN+"/{sample}_1.fastq.gz",
                r2=SAMPLE_DIR_ORIGIN+"/{sample}_2.fastq.gz"
        output:
                zipf="{result_path}/FastQC/{sample}_fastqc.zip",
                html="{result_path}/FastQC/{sample}_fastqc.html"
        params:
                "{result_path}/FastQC/}"
        shell:
                "fastqc -t 14 {input.r1} {input.r2} -o {params}"

rule TRIMMOMATIC:
        input:
                r1=SAMPLE_DIR_ORIGIN+"/{sample}_1.fastq.gz",
                r2=SAMPLE_DIR_ORIGIN+"/{sample}_2.fastq.gz"
        output:
                out1=SAMPLE_DIR+"/{sample}/{sample}_1.trimmed.fastq.gz",
                out2=SAMPLE_DIR+"/{sample}/{sample}_1.untrimmed.fastq.gz",
                out3=SAMPLE_DIR+"/{sample}/{sample}_2.trimmed.fastq.gz",
                out4=SAMPLE_DIR+"/{sample}/{sample}_2.untrimmed.fastq.gz"
                #out5=SAMPLE_DIR+"/{sample}/{sample}_single.trimmed.fastq.gz"
        params:
                pe_ref=PE_REF
        shell:
                "trimmomatic PE -threads 14 {input.r1} {input.r2} {output.out1} {output.out2} {output.out3} {output.out4} ILLUMINACLIP:{params.pe_ref}:2:30:10"

rule BwaMem:
	input:
		r1=SAMPLE_DIR+"/{sample}/{sample}_1.trimmed.fastq.gz",
		r2=SAMPLE_DIR+"/{sample}/{sample}_2.trimmed.fastq.gz"
	output:
		bam="{sample}_sorted.bam"
	params:
		bwa_t="10",
		bwa_K="100000000",
		bwa_v="3",
		bwa_Y=REF,
		bwa_R="'@RG\\tID:foo\\tLB:bar\\tPL:illumina\\tPU:illumina\\tSM:{sample}'"
	shell:
		"bwa mem -t {params.bwa_t} -K {params.bwa_K} -v {params.bwa_v} -Y {params.bwa_Y} "
		"{input.r1} {input.r2} -R {params.bwa_R} | "
		"samtools view -F 0x800 -1 -h -q 0 - | "
		"samtools sort -o {output.bam}"

rule Picard:
	input:
		rules.BwaMem.output.bam
	output:
		duplbam="{sample}_sorted.dupl.bam"
	params:
		met="{sample}_dupl.metrics.txt",
		VS="SILENT",
		ODPD="2500",
		ASO="coordinate",
		RD="true"
	shell:
		"picard MarkDuplicates -INPUT {input} -OUTPUT {output.duplbam} -METRICS_FILE {params.met} "
		"-VALIDATION_STRINGENCY {params.VS} -OPTICAL_DUPLICATE_PIXEL_DISTANCE {params.ODPD} "
		"-ASSUME_SORT_ORDER {params.ASO} -REMOVE_DUPLICATES {params.RD}"

rule BaseRecalibration:
	input:
		rules.Picard.output.duplbam
	output:
		recal="{sample}_recal.table"
	params:
		ref=REF,
		snp=KNOWN_VARIANTS,
		indel=KNOWN_INDELS
	shell:
		"gatk BaseRecalibrator -R {params.ref} --use-original-qualities --known-sites {params.snp} "
		"--known-sites {params.indel} -I {input} -O {output}"

rule ApplyBqsr:
	input:
		bqsr=rules.BaseRecalibration.output.recal,
		bam=rules.Picard.output.duplbam
	output:
		bam="{sample}_dupl.recalib.bam"
	params:
		ref=REF
	shell:
		"gatk ApplyBQSR --add-output-sam-program-record -R {params.ref} "
		"-bqsr {input.bqsr} --use-original-qualities --static-quantized-quals 10 "
		"--static-quantized-quals 20 --static-quantized-quals 30 "
		"-I {input.bam} -O {output.bam}"

rule Mutect2:
	input:
		normal=lambda wildcards: expand(f"/{config[wildcards.t_sample]['normal']}_dupl.recalib.bam",result_path=RESULT_DIR),
		tumor="{sample}_dupl.recalib.bam"
	output:
		vcf="{sample}_somatic.vcf",
		stats="{sample}_somatic.vcf.stats"
	params:
		ref=REF,
		normal=lambda wildcards: f"{config[wildcards.t_sample]['normal']}",
		germline=MUTECT_GERMLINE	
		#pon
	shell:
		"gatk Mutect2 -R {params.ref} -I {input.tumor} -I {input.normal} -normal {params.normal} -germline-resource {params.germline} -O {output.vcf}"

rule GetPileupSummaries:
	input:
		rules.ApplyBqsr.output.bam
	output:
		"{sample}_pileup.table"
	params:
		germline=MUTECT_GERMLINE,
		bed=BED
	shell:
		"gatk GetPileupSummaries -I {input} -V {params.germline} -L {params.bed} -O {output}"

rule CalculateConta:
	input:
		normal=lambda wildcards: expand(f"{config[wildcards.t_sample]['normal']}_pileup.table",result_path=RESULT_DIR),
		tumor="{t_sample}_pileup.table"
	output:
		conta="{t_sample}_contam.table",
		tumor_seg="{t_sample}_segments.table"
	shell:
		"gatk CalculateContamination -I {input.tumor} -matched {input.normal} -O {output.conta} --tumor-segmentation {output.tumor_seg}"

rule FilterMutectCalls:
        input:
                vcf=rules.Mutect2.output.vcf,
		stats=rules.Mutect2.output.stats,
                conta=rules.CalculateConta.output.conta,
		tumor_seg=rules.CalculateConta.output.tumor_seg
        output:
                vcf="{sample}_somatic.filtered.vcf"
	params:
		ref=REF,
		bed=BED
	shell:
		"gatk FilterMutectCalls -V {input.vcf} -R {params.ref} -L {params.bed} -O {output.vcf} "
		"--contamination-table {input.conta} --stats {input.stats} --tumor-segmentation {input.tumor_seg}"

rule PASSMUTATION:
	input:
		vcf=rules.FilterMutectCalls.output.vcf
	output:
		vcf="{t_sample}_somatic.filtered.PASS.vcf"
	shell:
		"bcftools view -f PASS {input.vcf} > {output}"

rule VCF2MAF:
        input:
                rules.PASSMUTATION.output.vcf
        output:
                "{t_sample}_somatic.filtered.PASS.maf"
	params:
		vcf2maf_tool="vcf2maf.pl",

		vep_tool="ensembl-vep",
		vep_cache="vep_cache",
		ref=REF,
		filtervcf=KNOWN_VARIANTS,
		tumor="{t_sample}",
		normal=lambda wildcards: f"{config[wildcards.t_sample]['normal']}"
	shell:
		#"perl {params.vcf2maf_tool} --input-vcf {input} --output-maf {output} "
		"{params.vcf2maf_tool} --input-vcf {input} --output-maf {output} "
		"--ref-fasta {params.ref} --vep-path {params.vep_tool} --vep-data {params.vep_cache} "
		"--tumor-id {params.tumor} --vcf-tumor-id {params.tumor} "
		"--normal-id {params.normal} --vcf-normal-id {params.normal} "
		"filter-vcf {params.filtervcf}"
rule CNVkit:
	input:
		normal=lambda wildcards: expand(f"{config[wildcards.t_sample]['normal']}_dupl.recalib.bam",result_path=RESULT_DIR),
		tumor="/{t_sample}_dupl.recalib.bam"
	output:
		cnn_file="{t_sample}_my_reference.cnn",
		cnn_dir="cnn_directory"
	params:
		cnvkit_tool="cnvkit.py",
		genebed="exome_regions.bed",
		anno="refFlat.txt",
		access="access.hg19.bed",
		ref=REF
	shell:
		"{params.cnvkit_tool} batch {input.tumor} --normal {input.normal} --fasta {params.ref} "
		"--access {params.access} --annotate {params.anno} "
		"--targets {params.genebed} --fasta {params.ref} "
		"--output-reference {output.cnn_file} --output-dir {output.cnn_dir} -p 40 --scatter --diagram"

