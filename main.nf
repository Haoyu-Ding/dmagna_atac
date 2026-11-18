nextflow.enable.dsl=2

/*
 * D. magna paired-end ATAC-seq pipeline.
 *
 * The sample sheet is the only source of sample selection. The validator
 * normalizes paths and permits repeated sample_id rows for lane merging.
 */

params.samplesheet = null
params.reference_fasta = null
params.bwa_index_dir = null
params.mitochondrial_contig = null
params.macs_genome_size = null
params.outdir = 'results'
params.container = 'ghcr.io/haoyu-ding/dmagna_atac:0.1.0'
params.min_length = 20
params.adapter_r1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
params.adapter_r2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
params.extra_cutadapt = ''
params.extra_bwa = ''
params.extra_macs = ''


process VALIDATE_INPUT {
    tag 'sample-sheet'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path samplesheet

    output:
    path 'samplesheet.used.csv', emit: validated_samplesheet
    path 'input_files.tsv', emit: input_files
    path 'input_validation_report.txt', emit: validation_report

    script:
    """
    python /opt/dmagna-atac/bin/validate_samplesheet.py \
        --input ${samplesheet} \
        --validated-output samplesheet.used.csv \
        --input-files-output input_files.tsv \
        --report-output input_validation_report.txt
    """
}


process RAW_FASTQC {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/qc/raw/${meta.sample_id}", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.sample_id}.raw_fastqc"), emit: raw_qc

    script:
    def read1_files = read1 instanceof List ? read1 : [read1]
    def read2_files = read2 instanceof List ? read2 : [read2]
    def read1_args = read1_files.collect { "'${it.name}'" }.join(' ')
    def read2_args = read2_files.collect { "'${it.name}'" }.join(' ')
    """
    mkdir -p ${meta.sample_id}.raw_fastqc
    fastqc --threads ${task.cpus} --outdir ${meta.sample_id}.raw_fastqc ${read1_args} ${read2_args}
    """
}


process TRIM_READS {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/qc/trimmed/${meta.sample_id}", mode: 'copy', pattern: '*.json'

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.sample_id}.trimmed.R1.fastq.gz"), path("${meta.sample_id}.trimmed.R2.fastq.gz"), emit: trimmed
    path "${meta.sample_id}.cutadapt.json", emit: cutadapt_report

    script:
    def read1_files = read1 instanceof List ? read1 : [read1]
    def read2_files = read2 instanceof List ? read2 : [read2]
    def read1_args = read1_files.collect { "'${it.name}'" }.join(' ')
    def read2_args = read2_files.collect { "'${it.name}'" }.join(' ')
    """
    set -euo pipefail
    cat ${read1_args} > ${meta.sample_id}.merged.R1.fastq.gz
    cat ${read2_args} > ${meta.sample_id}.merged.R2.fastq.gz

    cutadapt \
        --cores ${task.cpus} \
        --pair-filter=any \
        --minimum-length ${params.min_length} \
        -a ${params.adapter_r1} \
        -A ${params.adapter_r2} \
        --json=${meta.sample_id}.cutadapt.json \
        ${params.extra_cutadapt} \
        -o ${meta.sample_id}.trimmed.R1.fastq.gz \
        -p ${meta.sample_id}.trimmed.R2.fastq.gz \
        ${meta.sample_id}.merged.R1.fastq.gz \
        ${meta.sample_id}.merged.R2.fastq.gz
    """
}


process TRIMMED_FASTQC {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/qc/trimmed/${meta.sample_id}", mode: 'copy', pattern: '*_fastqc.*'

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.sample_id}.trimmed_fastqc"), emit: trimmed_qc

    script:
    """
    mkdir -p ${meta.sample_id}.trimmed_fastqc
    fastqc --threads ${task.cpus} --outdir ${meta.sample_id}.trimmed_fastqc ${read1} ${read2}
    """
}


process ALIGN_READS {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(read1), path(read2)
    path reference
    path bwa_index_files

    output:
    tuple val(meta), path("${meta.sample_id}.sorted.bam"), path("${meta.sample_id}.sorted.bam.bai"), emit: aligned

    script:
    def reference_name = reference.name
    """
    set -euo pipefail
    bwa-mem2 mem \
        -t ${task.cpus} \
        ${params.extra_bwa} \
        ${reference_name} \
        ${read1} \
        ${read2} \
        | samtools sort -@ ${task.cpus} -o ${meta.sample_id}.sorted.bam -
    samtools index -@ ${task.cpus} ${meta.sample_id}.sorted.bam
    """
}


process FILTER_BAM {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/filtered_bam", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(sorted_bam), path(sorted_bai)

    output:
    tuple val(meta), path("${meta.sample_id}.filtered.bam"), path("${meta.sample_id}.filtered.bam.bai"), emit: filtered
    path "${meta.sample_id}.filtered.flagstat.txt", emit: flagstat

    script:
    """
    set -euo pipefail
    samtools collate -@ ${task.cpus} -o ${meta.sample_id}.collated.bam ${sorted_bam}
    samtools fixmate -@ ${task.cpus} -m ${meta.sample_id}.collated.bam ${meta.sample_id}.fixmate.bam
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}.fixmate.sorted.bam ${meta.sample_id}.fixmate.bam
    samtools markdup -@ ${task.cpus} -r ${meta.sample_id}.fixmate.sorted.bam ${meta.sample_id}.deduplicated.bam

    samtools view -h -f 3 -F 3852 ${meta.sample_id}.deduplicated.bam \
        | awk -v mitochondrial_contig='${params.mitochondrial_contig}' \
            'BEGIN {OFS="\\t"} /^@/ {print; next} \$3 != mitochondrial_contig && \$7 != mitochondrial_contig {print}' \
        | samtools view -b -o ${meta.sample_id}.filtered.bam -
    samtools index -@ ${task.cpus} ${meta.sample_id}.filtered.bam
    samtools flagstat -@ ${task.cpus} ${meta.sample_id}.filtered.bam > ${meta.sample_id}.filtered.flagstat.txt
    """
}


process CALL_PEAKS {
    tag "${meta.sample_id}"

    publishDir "${params.outdir}/peaks", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(filtered_bam), path(filtered_bai)

    output:
    tuple val(meta), path("${meta.sample_id}.narrowPeak"), path("${meta.sample_id}.summits.bed"), emit: peaks
    path "${meta.sample_id}.peak_report.txt", emit: peak_report

    script:
    """
    set -euo pipefail
    mkdir -p macs3
    macs3 callpeak \
        -t ${filtered_bam} \
        -f BAMPE \
        -n ${meta.sample_id} \
        --outdir macs3 \
        --gsize ${params.macs_genome_size} \
        ${params.extra_macs}
    cp macs3/${meta.sample_id}_peaks.narrowPeak ${meta.sample_id}.narrowPeak
    cp macs3/${meta.sample_id}_summits.bed ${meta.sample_id}.summits.bed
    cp macs3/${meta.sample_id}_peaks.xls ${meta.sample_id}.peak_report.txt
    """
}


process MULTIQC_SUMMARY {
    tag 'run-level-qc'

    publishDir "${params.outdir}/qc/multiqc", mode: 'copy', pattern: '*'

    input:
    path raw_qc
    path trimmed_qc
    path cutadapt_reports
    path alignment_reports
    path peak_reports

    output:
    path 'multiqc', emit: multiqc

    script:
    """
    mkdir -p multiqc
    multiqc . \
        --force \
        --outdir multiqc \
        --filename multiqc_report.html
    """
}


process SOFTWARE_VERSIONS {
    tag 'software-versions'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    val trigger

    output:
    path 'software_versions.yml'

    script:
    """
    {
        echo 'fastqc: "'\$(fastqc --version 2>&1 | head -n 1)'"'
        echo 'cutadapt: "'\$(cutadapt --version 2>&1 | head -n 1)'"'
        echo 'bwa-mem2: "'\$(bwa-mem2 version 2>&1 | head -n 1)'"'
        echo 'samtools: "'\$(samtools --version 2>&1 | head -n 1)'"'
        echo 'macs3: "'\$(macs3 --version 2>&1 | head -n 1)'"'
        echo 'multiqc: "'\$(multiqc --version 2>&1 | head -n 1)'"'
    } > software_versions.yml
    """
}


workflow {
    if (!params.samplesheet) {
        error "Missing --samplesheet"
    }
    if (!params.reference_fasta) {
        error "Missing --reference_fasta: provide the fixed D. magna reference FASTA"
    }
    if (!params.bwa_index_dir) {
        error "Missing --bwa_index_dir: provide the directory containing the BWA-MEM2 index files"
    }
    if (!params.mitochondrial_contig) {
        error "Missing --mitochondrial_contig: provide the mitochondrial contig name for this reference"
    }
    if (!params.macs_genome_size) {
        error "Missing --macs_genome_size: provide the effective genome size selected for this reference"
    }

    samplesheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
    reference_ch = Channel.fromPath(params.reference_fasta, checkIfExists: true).first()
    bwa_index_ch = Channel.fromPath(
        "${params.bwa_index_dir}/*",
        type: 'file',
        checkIfExists: true,
    ).collect()

    VALIDATE_INPUT(samplesheet_ch)

    sample_rows = VALIDATE_INPUT.out.validated_samplesheet
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id: row.sample_id,
                condition: row.condition ?: '',
                replicate: row.replicate ?: '',
            ]
            tuple(meta, file(row.fastq_1), file(row.fastq_2))
        }
        .groupTuple(by: 0)

    RAW_FASTQC(sample_rows)
    TRIM_READS(sample_rows)
    TRIMMED_FASTQC(TRIM_READS.out.trimmed)
    ALIGN_READS(TRIM_READS.out.trimmed, reference_ch, bwa_index_ch)
    FILTER_BAM(ALIGN_READS.out.aligned)
    CALL_PEAKS(FILTER_BAM.out.filtered)

    MULTIQC_SUMMARY(
        RAW_FASTQC.out.raw_qc.map { meta, report_dir -> report_dir }.collect(),
        TRIMMED_FASTQC.out.trimmed_qc.map { meta, report_dir -> report_dir }.collect(),
        TRIM_READS.out.cutadapt_report.collect(),
        FILTER_BAM.out.flagstat.collect(),
        CALL_PEAKS.out.peak_report.collect(),
    )
    SOFTWARE_VERSIONS(Channel.value(true))
}
