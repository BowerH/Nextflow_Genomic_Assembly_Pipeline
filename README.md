# Nextflow_Genomic_Assembly_Pipeline

A bioinformatics workflow for processing sequencing data from SRA accession to assembly quality assessment.

**Made for Workflow Exercise in BIOL7210**

## Overview
**Workflow Steps**:
1. Prefetch SRA data (`sra-tools`)
2. Quality trimming (`fastp`)
3. Genome assembly (`SPAdes`)
4. Read filtering (`SeqKit`) *in parallel with assembly*
5. Assembly assessment (`QUAST`)

## Prerequisites
- [Nextflow](https://www.nextflow.io/) ≥22.10.8
- [Conda](https://docs.conda.io/en/latest/) ≥4.14
- 16GB+ RAM recommended

## Installation

```
git clone https://github.com/BowerH/Nextflow_Genomic_Assembly_Pipeline.git
cd Nextflow_Genomic_Assembly_Pipeline
```

## Usage
**Basic Execution**:
```
conda create -n nf -c bioconda nextflow -y
conda activate nf
nextflow run fullpipeline.nf -profile conda
```

Execution of mini dataset to ensure nextflow is working:

```
conda create -n nf -c bioconda nextflow -y
conda activate nf
nextflow run minipipeline.nf -profile conda
```

## Dependencies
```
Managed automatically via Conda:
| Tool       | Version  | Conda Channel  |
|------------|----------|----------------|
| sra-tools  | latest   | bioconda       |
| fastp      | 0.23.4   | bioconda       |
| SPAdes     | 3.15.5   | bioconda       |
| SeqKit     | 2.6.1    | bioconda       |
| QUAST      | 5.2.0    | bioconda       |

```

**Custom Parameters**:

```
nextflow run main.nf -profile conda
--srr_id "SRRXXXXXXX"
--outdir "custom_results"
```

## Output Structure

```
results/
├── quast/ # QUAST HTML reports
fastp_output/ # Trimmed reads & QC reports
spades_output/ # Assembly contigs
seqkit_raw_output/ # Filtered raw reads
data/ # Downloaded FASTQ files
```
