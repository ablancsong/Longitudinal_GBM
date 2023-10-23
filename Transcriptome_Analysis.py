
PATH="Path directory"

REF="ucsc.hg19.fasta"
PE_REF="TruSeq2-PE.fa"
STAR_DIR="StarIndex/"
GTF="gencode.v19.annotation.gtf"
HG19_REF_FLAT="refFlat.txt"

rule FASTQC:
	input:
		r1="{sample}_1.fastq.gz",
		r2="{sample}_2.fastq.gz"
	output:
		zipf="{sample}_fastqc.zip",
		html="{sample}_fastqc.html"
	shell:
		"fastqc -t 8 {input.r1} {input.r2} -o {output_directory}"

rule TRIMMOMATIC:
	input:
		r1="{sample}_1.fastq.gz",
		r2="{sample}_2.fastq.gz"
	output:
		out1="{sample}_1.trimmed.fastq.gz",
		out2="{sample}_1.untrimmed.fastq.gz",
		out3="{sample}_2.trimmed.fastq.gz",
		out4="{sample}_2.untrimmed.fastq.gz",
	params:
		pe_ref=PE_REF
	shell:
		"trimmomatic PE -threads 14 {input.r1} {input.r2} {output.out1} {output.out2} {output.out3} {output.out4} ILLUMINACLIP:{params.pe_ref}:2:30:10"

rule STAR:
	input:
		r1="{sample}_1.trimmed.fastq.gz",
		r2="{sample}_2.trimmed.fastq.gz"
	output:
		bam="{sample}_Aligned.sortedByCoord.out.bam",
		counts="{sample}_ReadsPerGene.out.tab"
	params:
		th="4",
		read_length="100",
		prefix="{sample}_",
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
		bai="{sample}_Aligned.sortedByCoord.out.bam.bai"
	shell:
	#	"samtools index -b -@ 14 {input.bam} -o {output.bai}"
		"samtools index -b {input.bam}"

rule htseq:
	input:
		bam=rules.STAR.output.bam,
		bai=rules.samtools_index.output.bai
	output:
		readcounts="{sample}_htseq.counts.out"
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
		standard="{sample}_stringtie.gtf",
		abundance="{sample}_stringtie.abund.gtf",
	params:
		th="8",
		gtf=GTF
	shell:
		"stringtie {input.bam} -p 8 -o {output.standard} -A {output.abundance} -G {params.gtf}"
rule STARFusion:
	input:
		r1="{sample}_1.trimmed.fastq.gz",
		r2="{sample}_2.trimmed.fastq.gz"
	params:
		genome="GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir",
		result_path="Result_path"
	shell:
		"STAR-Fusion --left_fq {input.r1} --right_fq {input.r2} "
		"--CPU 24 --output_dir {params.result_path} --genome_lib_dir {params.genome}"
