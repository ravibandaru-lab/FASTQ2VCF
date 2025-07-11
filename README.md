# <img alt="a bear" src="https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEg1ypWFO_NxVSMG2lbfE-gqHb8FNIF6nwjMG-A1GkTJJvTJTKhqe1fjmbQ_82O4SHPL_cqUBdT-vkcBXG1gjMC63bMHUv6wKCbKRY170aKXfukHlumOFg198kMJoEy7NJKVuEdiGITHmaYn/s800/animal_bear_character.png" height="60"> FASTQ2VCF

A Snakemake-based pipeline for variant calling that processes paired-end FASTQ files and generates filtered VCF files using standard bioinformatics tools (e.g., BWA, Samtools, GATK).

> [!CAUTION]
> This pipeline calls variants separately in each sample. Make sure this fits your usecase. If you want to call variants across multiple samples simultaneously (i.e. Group1 vs. Group2), reach out to me to modify the pipeline.

## Installation

#### Prerequsites
1. You need to have `conda` installed on your system. Learn more about conda [here](https://www.anaconda.com).

#### Downloading the Repository
```bash
$ git clone https://github.com/ravibandaru-lab/FASTQ2VCF
$ cd FASTQ2VCF
$ conda env create -f environment.yaml
$ conda activate FASTQ2VCF
```

#### Preparing the Pipeline

Create a folder-structure as follows:
```text
  root
	└── adapter
		└── TruSeq3-PE-2.fa
	└── fastq
		└── {SAMPLES}.R1.fastq.gz
		└── {SAMPLES}.R2.fastq.gz
	└── gene_annotation
		└── Gene_Annotation.gtf
	└── reference
		└── hg38.fa
	└── dag.pdf
	└── environment.yml
	└── LICENSE
	└── README
	└── Snakefile
```
When I was testing the pipeline:
- `Gene_Annotation.gtf`: I downloaded the comprehensive gene annotation, `gencode.v48.annotation.gtf`, from GENCODE [here](https://www.gencodegenes.org/human/).
- `hg38.fa`: I downloaded the reference genome for hg38 [here](https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/). Then `gunzip` to get .fa file.

#### Index the Genomes
```bash
cd reference
samtools faidx hg38.fa
bwa index hg38.fa
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict
cd ../
```

#### Modify the Snakefile
Go to the `Snakefile` and change `THREADS`, `GENE_ANNOTATION`, and `BWA_CORES` as per your constraints. I run on a `64` core machine, so my settings reflect that.

## Running the Pipeline
```bash
snakemake --cores {cores}
```

## Pipeline Diagram
See `dag.pdf`

## Pipeline Outputs
The FASTQ2VCF pipeline generates the following output files for each sample and reference genome combination:

- Raw FASTQ QC: `fastqc/pre/{sample}.R1_fastqc.html`, `fastqc/pre/{sample}.R2_fastqc.html`
- Post-trimming FASTQC: `fastqc/post/{sample}.R1_paired_fastqc.html`, `fastqc/post/{sample}.R2_paired_fastqc.html`

> [!IMPORTANT]
> FASTQ trimming was done with Trimmomatic using parameters `HEADCROP:6 ILLUMINACLIP:adapter/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36`. Make sure this fits your usecase.

- Paired and unpaired reads from Trimmomatic: `trimmed/{sample}.R1_paired.fastq.gz`, `trimmed/{sample}.R1_unpaired.fastq.gz`, `trimmed/{sample}.R2_paired.fastq.gz`, `trimmed/{sample}.R2_unpaired.fastq.gz`

> [!IMPORTANT]
> This pipeline currently supports paired-end sequencing data only. For each sample, both {sample}.R1.fastq.gz and {sample}.R2.fastq.gz files are required. Single-end data or unpaired reads will not be processed for alignment. Make sure this fits your usecase.

- Raw Aligned BAM: `bam/{sample}.{ref_genome}.raw.bam`
- Sorted and Indexed BAM: `bam/{sample}.{ref_genome}.sorted.bam`, `bam/{sample}.{ref_genome}.sorted.bam.bai`
- Filtered BAM (MAPQ ≥ 30; autosomes + chrX only): `bam/{sample}.{ref_genome}.filtered.bam`, `bam/{sample}.{ref_genome}.filtered.bam.bai`
> [!IMPORTANT]
> This filtered .bam is the one used for feature counting and variant calling. Make sure this fits your usecase.
- Raw Gene-Level Counts (via `featureCounts`): `feature_counts/{sample}.{ref_genome}.txt`
- `featureCounts` summary file: `feature_counts/{sample}.{ref_genome}.txt.summary`
- Custom Statistics: `feature_counts/{sample}.{ref_genome}.stats.txt`
> [!NOTE]
> Currently, the custom statistics reported are total genes, average length of gene, and average read count per gene.
- VCF file with SNP calls: `vcf/{sample}.{ref_genome}.calledSNPs.vcf`
> [!IMPORTANT]
> Previously `exactSNP` was used to call the SNPs, but was having a lot of trouble with that, so used `GATK HaplotypeCaller` instead. Only variants with a minimum Phred score of 30 (1 in 1000 chance of being wrong) are emitted. Make sure this fits your usecase.
## Contact
- Ravi Bandaru: ravi14.bandaru@gmail.com

## License
This project falls under an MIT license. See the included `LICENSE` file for details.
