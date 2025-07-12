# <img alt="a bear" src="https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEg1ypWFO_NxVSMG2lbfE-gqHb8FNIF6nwjMG-A1GkTJJvTJTKhqe1fjmbQ_82O4SHPL_cqUBdT-vkcBXG1gjMC63bMHUv6wKCbKRY170aKXfukHlumOFg198kMJoEy7NJKVuEdiGITHmaYn/s800/animal_bear_character.png" height="60"> FASTQ2VCF

A Snakemake-based pipeline for germline short variant discovery that processes paired-end FASTQ files and generates filtered VCF files using standard bioinformatics tools (e.g., BWA, Samtools, GATK).

## Installation

#### Prerequsites
1. You need to have `conda` installed on your system. Learn more about conda [here](https://www.anaconda.com).

#### Downloading the Repository
```bash
$ git clone https://github.com/ravibandaru-lab/FASTQ2VCF
$ cd FASTQ2VCF
$ conda env create -f environment.yaml
$ wget https://anaconda.org/bioconda/gatk4/4.6.2.0/download/noarch/gatk4-4.6.2.0-py310hdfd78af_0.tar.bz2
$ conda activate FASTQ2VCF
$ conda install gatk4-4.6.2.0-py310hdfd78af_0.tar.bz2
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
		└── known_sites.vcf.gz
		└── known_sites.vcf.gz.tbi
		└── known_indels.vcf.gz
		└── known_indels.vcf.gz.tbi
		└── 1000G.vcf.gz
		└── 1000G.vcf.gz.tbi
		└── hapmap.vcf.gz
		└── hapmap.vcf.gz.tbi
		└── omni.vcf.gz
		└── omni.vcf.gz.tbi
	└── dag.pdf
	└── environment.yml
	└── LICENSE
	└── README
	└── Snakefile
```
When I was testing the pipeline:
- `Gene_Annotation.gtf`: I downloaded the comprehensive gene annotation, `gencode.v48.annotation.gtf`, from GENCODE [here](https://www.gencodegenes.org/human/).
- `hg38.fa`: I downloaded the reference genome for hg38 [here](https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/). Then `gunzip` to get .fa file.
- `known_sites.vcf.gz(.tbi)` and `known_indels.vcf.gz(.tbi)`: 
```bash
# Download dbSNP (SNPs only)
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf -O reference/known_sites.vcf
bgzip reference/known_sites.vcf
tabix -p vcf reference/known_sites.vcf.gz

# Download Mills and 1000G gold indels
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O reference/known_indels.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi -O reference/known_indels.vcf.gz.tbi

#1000G
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz -O reference/1000G.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi -O reference/1000G.vcf.gz.tbi

# HapMap
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz -O reference/hapmap.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi -O reference/hapmap.vcf.gz.tbi

# Omni
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz -O reference/omni.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi -O reference/omni.vcf.gz.tbi
```

#### Modify the Snakefile
Go to the `Snakefile` and change `THREADS`, `GENE_ANNOTATION`, `BWA_CORES`, and `INTERVALS` as per your constraints. I run on a `64` core machine, so my settings reflect that.

## Running the Pipeline
```bash
snakemake --cores {cores}
```

## Pipeline Diagram
<img width="1082" height="1395" alt="Screenshot 2025-07-11 at 1 31 05 PM" src="https://github.com/user-attachments/assets/4009190f-e9a2-4eca-abb9-5a8b4bcedfe7" />


## Pipeline Outputs

1. **Quality Control**: FastQC analysis before and after trimming
2. **Read Preprocessing**: Trimmomatic adapter removal and quality trimming
3. **Alignment**: BWA-MEM alignment with read group information
4. **BAM Processing**: Duplicate marking and base quality score recalibration (BQSR)
5. **Variant Calling**: HaplotypeCaller in GVCF mode
6. **Joint Genotyping**: GenomicsDBImport followed by GenotypeGVCFs
7. **Variant Filtering**: VQSR for both SNPs and INDELs
8. **Gene Expression**: FeatureCounts for read quantification

The FASTQ2VCF pipeline generates the following output files for each sample and reference genome combination:

- Raw FASTQ QC: `fastqc/pre/{sample}.R1_fastqc.html`, `fastqc/pre/{sample}.R2_fastqc.html`
- Post-trimming FASTQC: `fastqc/post/{sample}.R1_paired_fastqc.html`, `fastqc/post/{sample}.R2_paired_fastqc.html`

> [!IMPORTANT]
> FASTQ trimming was done with Trimmomatic using parameters `HEADCROP:6 ILLUMINACLIP:adapter/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36`. Make sure this fits your usecase.

- Paired and unpaired reads from Trimmomatic: `trimmed/{sample}.R1_paired.fastq.gz`, `trimmed/{sample}.R1_unpaired.fastq.gz`, `trimmed/{sample}.R2_paired.fastq.gz`, `trimmed/{sample}.R2_unpaired.fastq.gz`

> [!IMPORTANT]
> This pipeline currently supports paired-end sequencing data only. For each sample, both {sample}.R1.fastq.gz and {sample}.R2.fastq.gz files are required. Single-end data or unpaired reads will not be processed for alignment. Make sure this fits your usecase.

- Raw Unaligned BAM: `ubam/{sample}.ubam.bam`
- Merged BAM: `alignment/{sample}.merged.bam`
- Deduplicated BAM: `dedup/{sample}.dedup.bam`, `dedup/{sample}.dedup.bam.bai`
- Recalibrated BAM (BQSR): `bqsr/{sample}.recalibrated.bam`, `bqsr/{sample}.recalibrated.bam.bai`
> [!IMPORTANT]
> This recalibrated .bam is the one used for feature counting and variant calling. 
- Raw Gene-Level Counts (via `featureCounts`): `feature_counts/{sample}.txt`
- `featureCounts` summary file: `feature_counts/{sample}.txt.summary`
- Custom Statistics: `feature_counts/{sample}.stats.txt`
> [!NOTE]
> Currently, the custom statistics reported are total genes, average length of gene, and average read count per gene. This is not normalized and does not account for GC content.
- GVCF Output for Each Sample: `gvcfs/{sample}.g.vcf.gz`, `gvcfs/{sample}.g.vcf.gz.tbi`
- Joint Calling Output: `joint_genotyping/genotyped.vcf.gz`, `joint_genotyping/genotyped.{chrom}.vcf.gz`
- VQSR Intermediate Files:
  - SNP-filtered VCF: `joint_genotyping/genotyped.filtered.snps.vcf.gz`
  - VQSR recalibration files: `joint_genotyping/recalibrate_SNP.*`, `joint_genotyping/recalibrate_INDEL.*`
- Final VQSR-Filtered VCF: `joint_genotyping/genotyped.filtered.vqsr.vcf.gz`
> [!IMPORTANT]
> Variant filtering uses VQSR (Variant Quality Score Recalibration) with standard GATK resources. SNPs are filtered at 99.5% sensitivity and INDELs at 99.0% sensitivity. The final output `genotyped.filtered.vqsr.vcf.gz` contains high-quality variants suitable for downstream analysis.
## Contact
- Ravi Bandaru: ravi14.bandaru@gmail.com

## License
This project falls under an MIT license. See the included `LICENSE` file for details.
