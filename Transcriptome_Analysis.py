import os,sys

#SAMPLE=os.listdir("/home/halim/HR_test/RNASeq/Fastq/SNU-E01T1-RNA")
SAMPLE=os.listdir("/home/halim/HR_test/RNASeq/Fastq/")
SAMPLE=[ '_'.join(i.split('_')[:1]) for i in SAMPLE ]
SAMPLE=list(set(SAMPLE))
# PATH
PATH="/home/halim/HR_test/RNASeq"
SAMPLE_DIR_ORIGIN=PATH+"/Fastq"
SAMPLE_DIR=PATH+"/Trim_Fastq"
RESULT_DIR=PATH+"/RESULT"
REF="/home/halim/HR_test/DB/ref_fasta/ucsc.hg19.fasta"
PE_REF="/home/halim/HR_test/DB/TruSeq2-PE.fa"
STAR_DIR="/home/halim/HR_test/DB/StarIndex/"
GTF="/home/halim/HR_test/DB/gencode.v19.annotation.gtf"
#GTF="/home/halim/HR_test/DB/Homo_sapiens.GRCh37.87.gtf"
HG19_REF_FLAT="/home/halim/HR_test/DB/refFlat.txt"

rule all:
        input:
                #expand("{result_path}/{sample}/{sample}_htseq.counts.out",result_path=RESULT_DIR,sample=SAMPLE),
                expand("{result_path}/{sample}/{sample}_stringtie.abund.gtf",result_path=RESULT_DIR,sample=SAMPLE)

rule FASTQC:
	input:
		r1=SAMPLE_DIR_ORIGIN+"/{sample}"+"/{sample}_1.fastq.gz",
		r2=SAMPLE_DIR_ORIGIN+"/{sample}"+"/{sample}_2.fastq.gz"
	output:
		zipf="{result_path}/FastQC/{sample}_fastqc.zip",
		html="{result_path}/FastQC/{sample}_fastqc.html"
	params:
		"{result_path}/FastQC/}"
	shell:
		"fastqc -t 8 {input.r1} {input.r2} -o {params}"

rule TRIMMOMATIC:
	input:
		r1=SAMPLE_DIR_ORIGIN+"/{sample}"+"/{sample}_1.fastq.gz",
		r2=SAMPLE_DIR_ORIGIN+"/{sample}"+"/{sample}_2.fastq.gz"
	output:
		out1=SAMPLE_DIR+"/{sample}/{sample}_1.trimmed.fastq.gz",
		out2=SAMPLE_DIR+"/{sample}/{sample}_1.untrimmed.fastq.gz",
		out3=SAMPLE_DIR+"/{sample}/{sample}_2.trimmed.fastq.gz",
		out4=SAMPLE_DIR+"/{sample}/{sample}_2.untrimmed.fastq.gz",
	params:
		pe_ref=PE_REF
	shell:
		"trimmomatic PE -threads 14 {input.r1} {input.r2} {output.out1} {output.out2} {output.out3} {output.out4} ILLUMINACLIP:{params.pe_ref}:2:30:10"

rule STAR:
	input:
		r1=SAMPLE_DIR+"/{sample}/{sample}_1.trimmed.fastq.gz",
		r2=SAMPLE_DIR+"/{sample}/{sample}_2.trimmed.fastq.gz"
	output:
		bam="{result_path}/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
		counts="{result_path}/{sample}/{sample}_ReadsPerGene.out.tab"
	params:
		th="4",
		read_length="100",
		prefix="{result_path}/{sample}/{sample}_",
		index_dir=STAR_DIR,
		ref=REF,
		gtf=GTF
	shell:
		"STAR --runThreadN 8 --genomeDir {params.index_dir} --readFilesIn {input.r1} {input.r2} "
		"--readFilesCommand zcat --outFileNamePrefix {params.prefix}  --outSAMtype BAM SortedByCoordinate "
		"--outBAMsortingThreadN 4 --quantMode GeneCounts --sjdbGTFfile {params.gtf}"

rule samtools_index:
	input:
		bam=rules.STAR.output.bam
	output:
		bai="{result_path}/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
	shell:
	#	"samtools index -b -@ 14 {input.bam} -o {output.bai}"
		"samtools index -b {input.bam}"

rule htseq:
	input:
		bam=rules.STAR.output.bam,
		bai=rules.samtools_index.output.bai
	output:
		readcounts="{result_path}/{sample}/{sample}_htseq.counts.out"
	params:
		gtf=GTF
	shell:
		"htseq-count -f bam -r pos -s no -t exon -i gene_name "
		"{input.bam} {params.gtf} > {output.readcounts}"

rule stringtie:
	input:
		bam=rules.STAR.output.bam,
		bai=rules.samtools_index.output.bai
	output:
		standard="{result_path}/{sample}/{sample}_stringtie.gtf",
		abundance="{result_path}/{sample}/{sample}_stringtie.abund.gtf",
	params:
		th="8",
		gtf=GTF
	shell:
		"stringtie {input.bam} -p 8 -o {output.standard} -A {output.abundance} -G {params.gtf}"
