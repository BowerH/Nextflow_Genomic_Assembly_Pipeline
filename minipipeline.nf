#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters with defaults
params.srr_id = 'SRR32426443'
params.outdir = 'results'

// Define the workflow
workflow {
    // Step 1: Prefetch SRR data
    prefetch_output = prefetch_data(params.srr_id)
    
    // Step 2: Run fastp for quality trimming
    fastp_output = fastp_trim(prefetch_output)
    
    // Step 3: Run SPAdes for assembly
    spades_output = spades_assemble(fastp_output)
    
    // Step 4: Filter assembly with seqkit
    // Path 2: Run seqkit on the raw reads (separate operation)
    seqkit_raw_output = seqkit_filter_raw(fastp_output)
    
    // Step 5: Run QUAST for assembly quality assessment
    quast_output = quast(spades_output)
    
    // View results
    quast_output.mix(seqkit_raw_output).view()
}

// Process: Prefetch SRR data
process prefetch_data {
    conda 'sra-tools'
    tag "${srr_id}"
    
    input:
    val srr_id
    
    output:
    tuple path("data/mini_1.fastq"), path("data/mini_2.fastq")
    
    script:
    """
    mkdir -p data/
    prefetch ${srr_id}
    fastq-dump --split-files --outdir data/ ${srr_id}
    # Subsample the FASTQ files (2000 lines = 500 reads for paired-end)
    #spades requires at least 10,000 reads for assembly
    head -n 40000 data/${srr_id}_1.fastq > data/mini_1.fastq
    head -n 40000 data/${srr_id}_2.fastq > data/mini_2.fastq
    """
}

// Process: Run fastp for quality trimming
process fastp_trim {
    conda 'bioconda::fastp=0.23.4'
    input:
    tuple path(r1_fastq), path(r2_fastq)
    
    output:
    tuple path("fastp_output/cleaned_R1.fastq"), path("fastp_output/cleaned_R2.fastq")
    
    script:
    """
    mkdir -p fastp_output/
    fastp -i ${r1_fastq} -I ${r2_fastq} \\
          -o fastp_output/cleaned_R1.fastq -O fastp_output/cleaned_R2.fastq \\
          --html fastp_output/fastp_report.html --json fastp_output/fastp_report.json
    """
}

// Process: Run SPAdes for assembly
process spades_assemble {
    conda 'bioconda::spades=3.15.5'
    input:
    tuple path(r1_clean), path(r2_clean)
    
    output:
    path "spades_output/contigs.fasta"
    
    script:
    """
    mkdir -p spades_output/
    spades.py --disable-gzip-output -o spades_output/ \\
              -1 ${r1_clean} -2 ${r2_clean} --memory 16 --phred-offset 33
    """
}

// Process: Filter assembly with seqkit
process seqkit_filter_raw {
    conda 'bioconda::seqkit=2.6.1'
    input:
    tuple path(r1_clean), path(r2_clean)
    
    output:
    path "seqkit_raw_output/filtered_reads.fasta"
    
    script:
    """
    mkdir -p seqkit_raw_output/
    seqkit fq2fa ${r1_clean} > seqkit_raw_output/r1.fasta
    seqkit fq2fa ${r2_clean} > seqkit_raw_output/r2.fasta
    cat seqkit_raw_output/r1.fasta seqkit_raw_output/r2.fasta > seqkit_raw_output/all_reads.fasta
    seqkit seq -m 100 seqkit_raw_output/all_reads.fasta > seqkit_raw_output/filtered_reads.fasta
    """
}

// Process: Run QUAST for assembly quality assessment
process quast {
    conda 'bioconda::quast=5.2.0'
    publishDir "${params.outdir}/quast", mode: 'copy'
    
    input:
    path filtered_assembly
    
    output:
    path "quast_output"
    
    script:
    """
    mkdir -p quast_output/
    quast.py ${filtered_assembly} -o quast_output/
    """
}