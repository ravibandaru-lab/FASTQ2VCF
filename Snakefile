import os
import glob
from tempfile import TemporaryDirectory

# Directories and files
READDIR = "fastq"
ADAPTERS = "adapter/TruSeq3-PE-2.fa"
GENE_ANNOTATION = "gene_annotation/gencode.v48.annotation.gtf"
KNOWN_SITES = [
    "reference/known_sites.vcf.gz",
    "reference/known_indels.vcf.gz"
]

PON = "reference/1000g_pon.vcf.gz"
GNOMAD_AF = "reference/gnomad_af_only.vcf.gz"

# Default settings
THREADS = 4
BWA_CORES = 8
REF_GENOME = "hg38"
FUNCOTATOR_DS = "reference/funcotator_dataSources.v1.8.hg38.20230908s"

# Reference genomes
REFERENCE_GENOMES = {
    "hg38": "reference/hg38.fa",
}

def get_ref_genome(ref=REF_GENOME):
    return REFERENCE_GENOMES.get(ref, f"reference/{ref}.fa")

# Samples and intervals
r1_files = glob.glob(f"{READDIR}/*.R1.fastq.gz")
SAMPLES = [os.path.basename(f).replace(".R1.fastq.gz", "") for f in r1_files]
INTERVALS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

rule all:
    input:
        expand("fastqc/pre/{sample}.R1_fastqc.html", sample=SAMPLES),
        expand("fastqc/pre/{sample}.R2_fastqc.html", sample=SAMPLES),
        expand("fastqc/post/{sample}.R1_paired_fastqc.html", sample=SAMPLES),
        expand("fastqc/post/{sample}.R2_paired_fastqc.html", sample=SAMPLES),
        expand("trimmed/{sample}.R1_paired.fastq.gz", sample=SAMPLES),
        expand("trimmed/{sample}.R2_paired.fastq.gz", sample=SAMPLES),
        expand("ubam/{sample}.ubam.bam", sample=SAMPLES),
        expand("alignment/{sample}.merged.bam", sample=SAMPLES),
        expand("dedup/{sample}.dedup.bam", sample=SAMPLES),
        expand("bqsr/{sample}.recalibrated.bam", sample=SAMPLES),
        expand("feature_counts/{sample}.txt", sample=SAMPLES),
        expand("feature_counts/{sample}.stats.txt", sample=SAMPLES),
        expand("mutect2/{sample}.unfiltered.vcf.gz", sample=SAMPLES),
        expand("mutect2/{sample}.filtered.vcf.gz", sample=SAMPLES),
        expand("mutect2/{sample}.funcotated.vcf.gz", sample=SAMPLES),
        expand("metrics/{sample}.alignment_summary_metrics.txt", sample=SAMPLES),
        "multiqc_report/multiqc_report.html"

rule bwa_index:
    input:
        fasta=get_ref_genome()
    output:
        amb=get_ref_genome().replace(".fa", ".fa.amb"),
        ann=get_ref_genome().replace(".fa", ".fa.ann"),
        bwt=get_ref_genome().replace(".fa", ".fa.bwt"),
        pac=get_ref_genome().replace(".fa", ".fa.pac"),
        sa=get_ref_genome().replace(".fa", ".fa.sa"),
        fai=get_ref_genome() + ".fai",
        dict=get_ref_genome().replace(".fa", ".dict")
    shell:
        """
        bwa index {input.fasta}
        samtools faidx {input.fasta}
        gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict}
        """

rule fastqc_pre:
    input:
        r1=f"{READDIR}/{{sample}}.R1.fastq.gz",
        r2=f"{READDIR}/{{sample}}.R2.fastq.gz"
    output:
        "fastqc/pre/{sample}.R1_fastqc.html",
        "fastqc/pre/{sample}.R1_fastqc.zip",
        "fastqc/pre/{sample}.R2_fastqc.html",
        "fastqc/pre/{sample}.R2_fastqc.zip"
    threads: THREADS
    shell:
        "mkdir -p fastqc/pre && fastqc -t {threads} -o fastqc/pre {input.r1} {input.r2}"

rule trimmomatic:
    input:
        f1=f"{READDIR}/{{sample}}.R1.fastq.gz",
        f2=f"{READDIR}/{{sample}}.R2.fastq.gz"
    output:
        "trimmed/{sample}.R1_paired.fastq.gz",
        "trimmed/{sample}.R1_unpaired.fastq.gz",
        "trimmed/{sample}.R2_paired.fastq.gz",
        "trimmed/{sample}.R2_unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        args=lambda wc: (
            f"HEADCROP:6 ILLUMINACLIP:{ADAPTERS}:2:30:10:2:keepBothReads "
            "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36"
        )
    threads: THREADS
    run:
        os.makedirs("trimmed", exist_ok=True)
        os.makedirs("logs/trimmomatic", exist_ok=True)

        with TemporaryDirectory() as tmp:
            shell(
                f"trimmomatic PE -threads {threads} "
                f"{input.f1} {input.f2} "
                f"{tmp}/R1P.fq.gz {tmp}/R1U.fq.gz {tmp}/R2P.fq.gz {tmp}/R2U.fq.gz "
                f"{params.args} 2>&1 | tee {log}"
            )
            for i, suffix in enumerate(["R1P", "R1U", "R2P", "R2U"]):
                shell(f"cp {tmp}/{suffix}.fq.gz {output[i]}")

rule fastqc_post:
    input:
        r1="trimmed/{sample}.R1_paired.fastq.gz",
        r2="trimmed/{sample}.R2_paired.fastq.gz"
    output:
        "fastqc/post/{sample}.R1_paired_fastqc.html",
        "fastqc/post/{sample}.R1_paired_fastqc.zip",
        "fastqc/post/{sample}.R2_paired_fastqc.html",
        "fastqc/post/{sample}.R2_paired_fastqc.zip"
    threads: THREADS
    shell:
        "mkdir -p fastqc/post && fastqc -t {threads} -o fastqc/post {input.r1} {input.r2}"

rule fastq_to_ubam:
    input:
        r1="trimmed/{sample}.R1_paired.fastq.gz",
        r2="trimmed/{sample}.R2_paired.fastq.gz"
    output:
        "ubam/{sample}.ubam.bam"
    log:
        "logs/fastq_to_ubam/{sample}.log"
    params:
        sample_name="{sample}",
        lib="lib1",
        platform="ILLUMINA",
        unit="unit1"
    shell:
        """
        mkdir -p ubam logs/fastq_to_ubam
        picard FastqToSam FASTQ={input.r1} FASTQ2={input.r2} OUTPUT={output} \
            READ_GROUP_NAME={params.sample_name} SAMPLE_NAME={params.sample_name} \
            LIBRARY_NAME={params.lib} PLATFORM_UNIT={params.unit} PLATFORM={params.platform} \
            2>&1 | tee {log}
        """

rule bwa_mem:
    input:
        fasta=get_ref_genome(),
        index=[get_ref_genome().replace(".fa", ".fa.amb"),
               get_ref_genome().replace(".fa", ".fa.ann"),
               get_ref_genome().replace(".fa", ".fa.bwt"),
               get_ref_genome().replace(".fa", ".fa.pac"),
               get_ref_genome().replace(".fa", ".fa.sa")],
        r1="trimmed/{sample}.R1_paired.fastq.gz",
        r2="trimmed/{sample}.R2_paired.fastq.gz"
    output:
        temp("alignment/{sample}.sam")
    threads: BWA_CORES
    shell:
        "mkdir -p alignment && bwa mem -t {threads} {input.fasta} {input.r1} {input.r2} > {output}"

rule merge_bam:
    input:
        fasta=get_ref_genome(),
        sam="alignment/{sample}.sam",
        ubam="ubam/{sample}.ubam.bam"
    output:
        "alignment/{sample}.merged.bam"
    log:
        "logs/merge_bam/{sample}.log"
    shell:
        """
        mkdir -p logs/merge_bam
        picard MergeBamAlignment R={input.fasta} UNMAPPED_BAM={input.ubam} ALIGNED_BAM={input.sam} \
            O={output} CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false \
            CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            ATTRIBUTES_TO_RETAIN=XS 2>&1 | tee {log}
        """

rule mark_duplicates:
    input:
        bam="alignment/{sample}.merged.bam"
    output:
        bam="dedup/{sample}.dedup.bam",
        bai="dedup/{sample}.dedup.bam.bai"
    log:
        "logs/mark_duplicates/{sample}.log"
    shell:
        """
        mkdir -p dedup logs/mark_duplicates
        gatk MarkDuplicatesSpark -I {input.bam} -O {output.bam} --TMP_DIR tmp 2>&1 | tee {log}
        samtools index {output.bam}
        """

rule base_recalibrator:
    input:
        bam="dedup/{sample}.dedup.bam",
        bai="dedup/{sample}.dedup.bam.bai",
        fasta=get_ref_genome(),
        known_sites=KNOWN_SITES
    output:
        "bqsr/{sample}.recal_data.table"
    log:
        "logs/bqsr/{sample}.base_recal.log"
    shell:
        """
        mkdir -p bqsr logs/bqsr
        gatk BaseRecalibrator -I {input.bam} -R {input.fasta} \
            --known-sites {input.known_sites[0]} --known-sites {input.known_sites[1]} \
            -O {output} 2>&1 | tee {log}
        """

rule apply_bqsr:
    input:
        bam="dedup/{sample}.dedup.bam",
        table="bqsr/{sample}.recal_data.table",
        fasta=get_ref_genome()
    output:
        bam="bqsr/{sample}.recalibrated.bam",
        bai="bqsr/{sample}.recalibrated.bam.bai"
    log:
        "logs/bqsr/{sample}.apply_bqsr.log"
    shell:
        """
        gatk ApplyBQSR -R {input.fasta} -I {input.bam} \
            --bqsr-recal-file {input.table} -O {output.bam} 2>&1 | tee {log}
        
        samtools index -b {output.bam} {output.bai}
        """

rule featurecounts:
    input:
        bam = "bqsr/{sample}.recalibrated.bam",
        bai = "bqsr/{sample}.recalibrated.bam.bai"
    output:
        counts = "feature_counts/{sample}.txt",
        summary = "feature_counts/{sample}.txt.summary"
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
        counts = "feature_counts/{sample}.txt"
    output:
        stats = "feature_counts/{sample}.stats.txt"
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

rule mutect2:
    input:
        bam="bqsr/{sample}.recalibrated.bam",
        bai="bqsr/{sample}.recalibrated.bam.bai",
        fasta=get_ref_genome(),
        pon=PON,
        gnomad=GNOMAD_AF
    output:
        vcf="mutect2/{sample}.unfiltered.vcf.gz",
        vcf_index="mutect2/{sample}.unfiltered.vcf.gz.tbi",
        f1r2="mutect2/{sample}.f1r2.tar.gz"
    log:
        "logs/mutect2/{sample}.log"
    shell:
        """
        mkdir -p mutect2 logs/mutect2
        gatk Mutect2 \
            -R {input.fasta} \
            -I {input.bam} \
            --panel-of-normals {input.pon} \
            --germline-resource {input.gnomad} \
            -O {output.vcf} \
            --f1r2-tar-gz {output.f1r2} \
            2>&1 | tee {log}
        """

rule learn_orientation_model:
    input:
        f1r2="mutect2/{sample}.f1r2.tar.gz"
    output:
        tar="mutect2/{sample}.read-orientation-model.tar.gz"
    log:
        "logs/mutect2/{sample}.orientation.log"
    shell:
        """
        gatk LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {output.tar} \
            2>&1 | tee {log}
        """

rule get_pileup_summaries:
    input:
        bam="bqsr/{sample}.recalibrated.bam",
        bai="bqsr/{sample}.recalibrated.bam.bai",
        ref="reference/hg38.fa",
        gnomad="reference/gnomad_af_only.vcf.gz"
    output:
        "mutect2/{sample}.pileups.table"
    params:
        intervals=" -L ".join(INTERVALS)
    log:
        "logs/mutect2/{sample}.pileups.log"
    shell:
        """
        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.gnomad} \
            -R {input.ref} \
            -L {params.intervals} \
            -O {output} \
            2>&1 | tee {log}
        """

rule calculate_contamination:
    input:
        pileups="mutect2/{sample}.pileups.table"
    output:
        contamination="mutect2/{sample}.contamination.table",
        segments="mutect2/{sample}.segments.table"
    log:
        "logs/mutect2/{sample}.contamination.log"
    shell:
        """
        gatk CalculateContamination \
            -I {input.pileups} \
            -O {output.contamination} \
            --tumor-segmentation {output.segments} \
            2>&1 | tee {log}
        """

rule filter_mutect_calls:
    input:
        vcf="mutect2/{sample}.unfiltered.vcf.gz",
        vcf_index="mutect2/{sample}.unfiltered.vcf.gz.tbi",
        fasta=get_ref_genome(),
        contamination="mutect2/{sample}.contamination.table",
        orientation="mutect2/{sample}.read-orientation-model.tar.gz"
    output:
        vcf="mutect2/{sample}.filtered.vcf.gz",
        vcf_index="mutect2/{sample}.filtered.vcf.gz.tbi"
    log:
        "logs/mutect2/{sample}.filter.log"
    shell:
        """
        gatk FilterMutectCalls \
            -V {input.vcf} \
            -R {input.fasta} \
            --contamination-table {input.contamination} \
            --ob-priors {input.orientation} \
            -O {output.vcf} \
            2>&1 | tee {log}
        """

rule funcotator:
    input:
        vcf="mutect2/{sample}.filtered.vcf.gz",
        fasta=get_ref_genome()
    output:
        "mutect2/{sample}.funcotated.vcf.gz"
    log:
        "logs/mutect2/{sample}.funcotator.log"
    shell:
        """
        gatk Funcotator \
            -R {input.fasta} \
            -V {input.vcf} \
            -O {output} \
            --output-file-format VCF \
            --data-sources-path {FUNCOTATOR_DS} \
            --ref-version hg38 \
            2>&1 | tee {log}
        """

rule alignment_metrics:
    input:
        bam="bqsr/{sample}.recalibrated.bam",
        bai="bqsr/{sample}.recalibrated.bam.bai",
        fasta=get_ref_genome()
    output:
        metrics="metrics/{sample}.alignment_summary_metrics.txt"
    log:
        "logs/picard_metrics/{sample}.alignment_metrics.log"
    shell:
        """
        mkdir -p metrics logs/picard_metrics
        picard CollectAlignmentSummaryMetrics \
            R={input.fasta} \
            I={input.bam} \
            O={output.metrics} \
            2>&1 | tee {log}
        """

rule multiqc:
    input:
        fastqc_pre=expand("fastqc/pre/{sample}.R1_fastqc.html", sample=SAMPLES),
        fastqc_post=expand("fastqc/post/{sample}.R1_paired_fastqc.html", sample=SAMPLES),
        alignment_metrics=expand("metrics/{sample}.alignment_summary_metrics.txt", sample=SAMPLES),
        feature_counts=expand("feature_counts/{sample}.txt.summary", sample=SAMPLES)
    output:
        html="multiqc_report/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        """
        mkdir -p multiqc_report logs/multiqc
        multiqc fastqc/pre fastqc/post metrics feature_counts -o multiqc_report \
            2>&1 | tee {log}
        """