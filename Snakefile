import os
import glob
import shutil
from tempfile import TemporaryDirectory

READDIR = "fastq"
THREADS = 4
ADAPTERS = "adapter/TruSeq3-PE-2.fa" #from https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE-2.fa
GENE_ANNOTATION = "gene_annotation/gencode.v48.annotation.gtf"
BWA_CORES = 8

r1_files = glob.glob(f"{READDIR}/*.R1.fastq.gz")
SAMPLES = [os.path.basename(f).replace(".R1.fastq.gz", "") for f in r1_files]

REF_GENOMES = ["hg38"]

def get_ref_genome(ref_genome):
    """Return path to reference genome"""
    ref_paths = {
        "hg38": "reference/hg38.fa",
    }
    return ref_paths.get(ref_genome, f"reference/{ref_genome}.fa")

rule all:
    input:
        # Pre-trim FastQC
        expand("fastqc/pre/{sample}.R1_fastqc.html", sample=SAMPLES),
        expand("fastqc/pre/{sample}.R2_fastqc.html", sample=SAMPLES),

        # Trimmed FASTQ files
        expand("trimmed/{sample}.R1_paired.fastq.gz", sample=SAMPLES),
        expand("trimmed/{sample}.R1_unpaired.fastq.gz", sample=SAMPLES),
        

        # Post-trim FastQC
        expand("fastqc/post/{sample}.R1_paired_fastqc.html", sample=SAMPLES),
        expand("fastqc/post/{sample}.R2_paired_fastqc.html", sample=SAMPLES),

        # BAM files
        expand("bam/{sample}.{ref_genome}.raw.bam", sample=SAMPLES, ref_genome=REF_GENOMES),
        expand("bam/{sample}.{ref_genome}.sorted.bam", sample=SAMPLES, ref_genome=REF_GENOMES),
        expand("bam/{sample}.{ref_genome}.sorted.bam.bai", sample=SAMPLES, ref_genome=REF_GENOMES),
        expand("bam/{sample}.{ref_genome}.filtered.bam", sample=SAMPLES, ref_genome=REF_GENOMES),
        expand("bam/{sample}.{ref_genome}.filtered.bam.bai", sample=SAMPLES, ref_genome=REF_GENOMES),

        # featureCounts raw count files
        expand("feature_counts/{sample}.{ref_genome}.txt", sample=SAMPLES, ref_genome=REF_GENOMES),

        # featureCounts summarized output (simplified + stats)
        expand("feature_counts/{sample}.{ref_genome}.stats.txt", sample=SAMPLES, ref_genome=REF_GENOMES),
        
        # VCF files from exactSNP
        expand("vcf/{sample}.{ref_genome}.calledSNPs.vcf", sample=SAMPLES, ref_genome=REF_GENOMES),

rule fastqc_pre:
    input:
        r1 = f"{READDIR}/{{sample}}.R1.fastq.gz",
        r2 = f"{READDIR}/{{sample}}.R2.fastq.gz"
    output:
        r1_html = "fastqc/pre/{sample}.R1_fastqc.html",
        r1_zip = "fastqc/pre/{sample}.R1_fastqc.zip",
        r2_html = "fastqc/pre/{sample}.R2_fastqc.html",
        r2_zip = "fastqc/pre/{sample}.R2_fastqc.zip"
    threads: THREADS
    shell:
        """
        mkdir -p fastqc/pre
        fastqc -t {threads} -o fastqc/pre {input.r1} {input.r2}
        """

rule trimmomatic:
    input:
        fastq = [f"{READDIR}/{{sample}}.R1.fastq.gz", f"{READDIR}/{{sample}}.R2.fastq.gz"]
    output:
        trim = [
            "trimmed/{sample}.R1_paired.fastq.gz",
            "trimmed/{sample}.R1_unpaired.fastq.gz",
            "trimmed/{sample}.R2_paired.fastq.gz",
            "trimmed/{sample}.R2_unpaired.fastq.gz",
        ],
        log = "logs/trimmomatic/{sample}.log"
    params:
        trimmomatic_args = lambda wildcards: ("HEADCROP:6 ILLUMINACLIP:" + ADAPTERS + ":2:30:10:2:keepBothReads "
                                            "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36")
    threads: THREADS
    run:
        import os
        os.makedirs("trimmed", exist_ok=True)
        os.makedirs("logs/trimmomatic", exist_ok=True)
        trimmomatic_args = params.trimmomatic_args
        sample = wildcards.sample
        
        with TemporaryDirectory() as tempdir:
            cmd = (
                f"trimmomatic PE -threads {threads} "
                f"{input.fastq[0]} {input.fastq[1]} "
                f"{tempdir}/R1P.fq.gz {tempdir}/R1U.fq.gz "
                f"{tempdir}/R2P.fq.gz {tempdir}/R2U.fq.gz "
                f"{trimmomatic_args} "
                f"2>&1 | tee {output.log}"
            )
            shell(cmd)
            
            for tmp, final in zip(
                [f"{tempdir}/R1P.fq.gz", f"{tempdir}/R1U.fq.gz", f"{tempdir}/R2P.fq.gz", f"{tempdir}/R2U.fq.gz"],
                output.trim
            ):
                shutil.copy(tmp, final)
                os.remove(tmp)

rule fastqc_post:
    input:
        r1 = "trimmed/{sample}.R1_paired.fastq.gz",
        r2 = "trimmed/{sample}.R2_paired.fastq.gz"
    output:
        r1_html = "fastqc/post/{sample}.R1_paired_fastqc.html",
        r1_zip = "fastqc/post/{sample}.R1_paired_fastqc.zip",
        r2_html = "fastqc/post/{sample}.R2_paired_fastqc.html",
        r2_zip = "fastqc/post/{sample}.R2_paired_fastqc.zip"
    threads: THREADS
    shell:
        """
        mkdir -p fastqc/post
        fastqc -t {threads} -o fastqc/post {input.r1} {input.r2}
        """

rule bwa_raw:
    input: 
        ref = lambda wildcards: get_ref_genome(wildcards.ref_genome),
        fastq = ["trimmed/{sample}.R1_paired.fastq.gz", "trimmed/{sample}.R2_paired.fastq.gz"]
    output: 
        bam = "bam/{sample}.{ref_genome,[^\\.]+}.raw.bam"
    log:
        log = "logs/bwa/{sample}.{ref_genome}.samblaster.log"
    params:
        extra = lambda wildcards: f"-R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA'"
    threads: BWA_CORES
    shell:
        """
        mkdir -p bam logs/bwa
        
        bwa mem -t {threads} {params.extra} {input.ref} {input.fastq} | \\
        samblaster 2> >(tee {log.log} >&2) | \\
        samtools view -Sb - > {output.bam}
        """

rule bam_sort:
    input:
        bam = "bam/{sample}.{ref_genome}.raw.bam"
    output:
        sorted_bam = "bam/{sample}.{ref_genome}.sorted.bam",
        index = "bam/{sample}.{ref_genome}.sorted.bam.bai"
    threads: THREADS
    shell:
        """
        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam}
        samtools index {output.sorted_bam}
        """

rule bam_filter:
    input:
        bam = "bam/{sample}.{ref_genome}.sorted.bam",
        bai = "bam/{sample}.{ref_genome}.sorted.bam.bai"
    output:
        filtered_bam = "bam/{sample}.{ref_genome}.filtered.bam",
        filtered_index = "bam/{sample}.{ref_genome}.filtered.bam.bai"
    threads: THREADS
    shell:
        """
        samtools view -@ {threads} -b -q 30 {input.bam} \
            chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
            chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
            chr21 chr22 chrX > {output.filtered_bam}
        samtools index {output.filtered_bam}
        """

rule featurecounts:
    input:
        bam = "bam/{sample}.{ref_genome}.filtered.bam",
        bai = "bam/{sample}.{ref_genome}.filtered.bam.bai"
    output:
        counts = "feature_counts/{sample}.{ref_genome}.txt",
        summary = "feature_counts/{sample}.{ref_genome}.txt.summary"
    params:
        gtf = GENE_ANNOTATION
    threads: THREADS
    shell:
        """
        mkdir -p feature_counts
        featureCounts \
            -T {threads} \
            -p -B -C \
            -t exon \
            -g gene_id \
            -a {params.gtf} \
            -o {output.counts} \
            {input.bam}
        """

rule summarize_featurecounts:
    input:
        counts = "feature_counts/{sample}.{ref_genome}.txt"
    output:
        stats = "feature_counts/{sample}.{ref_genome}.stats.txt"
    shell:
        """
        awk 'NR > 2 && $7 > 0 {{
            count++;
            length_sum += $6;
            count_sum += $7
        }}
        END {{
            if (count > 0) {{
                print "Total expressed genes:", count;
                print "Average gene length:", length_sum/count;
                print "Average count per gene:", count_sum/count
            }} else {{
                print "No expressed genes found"
            }}
        }}' {input.counts} > {output.stats}
        """

rule bam_to_vcf:
    input:
        bam = "bam/{sample}.{ref_genome}.filtered.bam",
        bai = "bam/{sample}.{ref_genome}.filtered.bam.bai",
        ref = lambda wildcards: get_ref_genome(wildcards.ref_genome)
    output:
        vcf = "vcf/{sample}.{ref_genome}.calledSNPs.vcf"
    threads: THREADS
    shell:
        """
        mkdir -p vcf {params.tmp_dir}

		gatk HaplotypeCaller \
			-R {input.ref} \
			-I {input.bam} \
			-O {output.vcf} \
			-stand-call-conf 30.0 \
			-ERC NONE
        """
