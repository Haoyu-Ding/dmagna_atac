# D. magna ATAC-seq Pipeline Design

## 1. Purpose

This document defines the design of the first version of an independent ATAC-seq analysis pipeline for *Daphnia magna*.

The pipeline will be developed with Nextflow, executed on the KU Leuven VSC, and run with version-controlled containerised software. The workflow is inspired by the general structure of `nf-core/atacseq`, but it is an independent implementation with a smaller and explicitly defined scope.

The first version will process an arbitrary number of paired-end ATAC-seq samples and produce quality-controlled alignment files and per-sample peak calls.

## 2. Scope of version 1

### Included

- Paired-end Illumina FASTQ input.
- A user-provided CSV sample sheet as the only supported input mode.
- An arbitrary number of samples.
- Input validation before computational analysis.
- Raw-read quality control.
- Adapter and quality trimming.
- Alignment to a fixed *D. magna* reference genome.
- BAM sorting and indexing.
- Filtering of unsuitable alignments.
- PCR duplicate and mitochondrial-read handling.
- Per-sample peak calling in paired-end mode.
- Per-sample and run-level QC summaries.
- Nextflow execution reports and a record of the validated input.

### Not included

- Automatic discovery of FASTQ files from a directory.
- Automatic interpretation of sample names or sequencing indexes.
- Automatic assignment of Fish, Control, treatment, or replicate status.
- Analysis of `Undetermined` reads.
- Use of I1/I2 index-read FASTQ files for peak calling.
- Differential accessibility analysis.
- Consensus peak construction.
- Peak annotation, motif analysis, and functional enrichment.
- Biological conclusions about Fish versus Control samples.

## 3. Input design principle

The workflow must not contain a fixed list of sample names or assume that an experiment contains four samples. The number of samples is determined by the number of valid rows in the sample sheet.

The sample sheet is the authoritative allow-list for a run:

- A file is analysed only if it appears in the sample sheet.
- Files not listed in the sample sheet are ignored.
- Sample identity and biological metadata come from the sample sheet, not from guesses based on filenames.
- R1/R2 pairing is checked during validation.

This design prevents unrelated samples, `Undetermined` reads, or additional demultiplexed libraries in the same directory from being analysed accidentally.

## 4. Sample sheet specification

The first version will accept a comma-separated file with the following columns:

```text
sample_id,fastq_1,fastq_2,condition,replicate
```

### Required columns

| Column | Meaning | Requirement |
|---|---|---|
| `sample_id` | Stable identifier for one biological sample or analysed library | Required; unique after lane/run grouping |
| `fastq_1` | Path to the R1 FASTQ file | Required; must exist and be gzipped |
| `fastq_2` | Path to the R2 FASTQ file | Required; must exist and be gzipped |

### Optional metadata columns

| Column | Meaning | Version 1 use |
|---|---|---|
| `condition` | Experimental condition, for example `fish` or `control` | Recorded and propagated as metadata; not used for peak calling decisions |
| `replicate` | Biological replicate identifier | Recorded and propagated as metadata; not used for statistical testing |

Additional columns may be added later, but unknown columns should generate a warning rather than silently changing workflow behaviour.

## 5. Sample, library, and lane semantics

In version 1, one `sample_id` represents one unit for which a final BAM file and a final peak file will be produced.

If the same biological sample was sequenced over multiple lanes or split into multiple FASTQ files, the rows must use the same `sample_id`:

```text
sample_A,/data/sample_A_L001_R1.fastq.gz,/data/sample_A_L001_R2.fastq.gz,fish,1
sample_A,/data/sample_A_L002_R1.fastq.gz,/data/sample_A_L002_R2.fastq.gz,fish,1
```

The pipeline will validate that the metadata for repeated `sample_id` values is consistent, then combine the corresponding read pairs before downstream sample-level analysis.

Different biological replicates must have different `sample_id` values, even when they belong to the same condition:

```text
fish_rep1,...,fish,1
fish_rep2,...,fish,2
```

The current project-specific assignment of `condition` and `replicate` must be confirmed from the experimental record before the final sample sheet is committed.

## 6. Example sample sheet for the current dataset

The following is an example only. The exact file paths will be filled in on the VSC after the data transfer.

```csv
sample_id,fastq_1,fastq_2,condition,replicate
fish_r4,/path/to/GC187803_CCTGTCAT-ATGACAGG_S17_L001_R1_001.fastq.gz,/path/to/GC187803_CCTGTCAT-ATGACAGG_S17_L001_R2_001.fastq.gz,fish,1
fish_r5,/path/to/GC187803_TTAGGTCG-CGACCTAA_S22_L001_R1_001.fastq.gz,/path/to/GC187803_TTAGGTCG-CGACCTAA_S22_L001_R2_001.fastq.gz,fish,2
control_r4,/path/to/GC187803_GTGCCATA-TATGGCAC_S20_L001_R1_001.fastq.gz,/path/to/GC187803_GTGCCATA-TATGGCAC_S20_L001_R2_001.fastq.gz,control,1
control_r1,/path/to/GC187803_GGCGTTAT-ATAACGCC_S19_L001_R1_001.fastq.gz,/path/to/GC187803_GGCGTTAT-ATAACGCC_S19_L001_R2_001.fastq.gz,control,2
```

The workflow itself will not contain these four names. They are only one possible sample sheet supplied by the user for one run.

The following files are intentionally not included:

- `CGTACGAA-TTCGTACG` sample files;
- `TGGATCAC-GTGATCCA` sample files;
- all `Undetermined` files;
- all `I1` and `I2` index-read files.

## 7. Input validation

Input validation is a separate first stage and must complete before FastQC or alignment starts.

The validation stage must check:

1. The sample sheet exists and has the required header columns.
2. Every non-empty `sample_id` is valid for use in output filenames.
3. Every `fastq_1` and `fastq_2` path exists.
4. FASTQ files use `.fastq.gz` or `.fq.gz` extensions.
5. Every sample has both R1 and R2.
6. R1 and R2 filenames belong to the same sequencing unit.
7. Repeated `sample_id` rows have consistent condition and replicate metadata.
8. No sample has conflicting or duplicated read-pair entries.
9. At least one valid sample remains after validation.

The pipeline should stop with a clear error for invalid input. It must not silently remove files or repair ambiguous pairings.

The validation implementation should use a sample-sheet schema, preferably through the Nextflow schema ecosystem, so required fields, file paths, and field formats are checked before the workflow channel is created.

## 8. Nextflow data model

After validation, the sample sheet is converted into a dynamic channel. Each emitted item contains metadata and a paired set of reads:

```text
meta, read1, read2
```

Conceptually, `meta` contains:

```text
sample_id
condition
replicate
```

All processes receive this dynamic tuple. No process should refer directly to a fixed sample name or a fixed number of samples.

The workflow therefore scales from one sample to many samples without changing the Nextflow code:

```text
validated samplesheet
        ↓
sample-level metadata channel
        ↓
FastQC → trimming → alignment → BAM processing → peak calling
```

Independent samples should run in parallel whenever VSC resources permit. The downstream aggregation steps collect results from the channel rather than assuming a fixed sample count.

## 9. Computational workflow

```mermaid
flowchart LR
    A["CSV sample sheet"] --> B["Validate files and metadata"]
    B --> C["Group repeated rows by sample_id"]
    C --> D["Raw-read QC"]
    D --> E["Adapter and quality trimming"]
    E --> F["Paired-end alignment"]
    F --> G["Sort and index BAM"]
    G --> H["Filter alignments"]
    H --> I["Remove duplicates and mitochondrial reads"]
    I --> J["MACS3 paired-end peak calling"]
    J --> K["Per-sample peaks and QC"]
    K --> L["Run-level summary"]
```

The first version will call peaks independently for each `sample_id`. It will not combine Fish or Control samples into consensus peaks until the biological replicate structure and the intended downstream analysis have been confirmed.

## 10. Expected outputs

The result directory should contain a sample-level structure similar to:

```text
results/
├── qc/
│   ├── raw/
│   ├── trimmed/
│   └── multiqc/
├── alignment/
│   ├── <sample_id>.sorted.bam
│   └── <sample_id>.sorted.bam.bai
├── filtered_bam/
│   ├── <sample_id>.filtered.bam
│   └── <sample_id>.filtered.bam.bai
├── peaks/
│   ├── <sample_id>.narrowPeak
│   ├── <sample_id>.summits.bed
│   └── <sample_id>.peak_report.txt
└── pipeline_info/
    ├── samplesheet.used.csv
    ├── input_files.tsv
    ├── input_validation_report.txt
    ├── software_versions.yml
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

The copied `samplesheet.used.csv` is important: it records the exact input selection for that run, including the files that were intentionally excluded from the directory.

## 11. Container and execution design

The workflow code and Dockerfile will be stored in GitHub. The pipeline will not run on the local computer.

The local workspace is used to:

- edit workflow code;
- edit the Dockerfile and configuration;
- record metadata and documentation;
- commit and push changes to GitHub.

The VSC is used to:

- clone or pull the repository;
- access the FASTQ and reference files;
- run Nextflow through the VSC scheduler;
- run the containerised tools;
- inspect and archive results.

On VSC, the container execution layer will be designed for Apptainer while the Dockerfile remains the software-environment definition. VSC documents Apptainer as the HPC-compatible runtime for Docker/OCI-derived images. See the [VSC container documentation](https://docs.vscentrum.be/compute/software/installing_software/containers.html).

## 12. First-version acceptance criteria

Version 1 is considered functional when it can:

- accept a sample sheet with an arbitrary number of valid paired-end samples;
- reject missing, unmatched, or ambiguous input files before analysis;
- process each valid sample without hard-coded sample names;
- produce one final filtered BAM and one peak result per sample;
- produce a run-level QC summary;
- record the exact validated input sample sheet;
- run reproducibly on the VSC with the pinned container environment;
- process the current four-sample dataset using a user-provided sample sheet.

## 13. Design decision

The project will use a samplesheet-first input model. Automatic FASTQ directory scanning will not be part of the first version. This keeps sample selection explicit, prevents accidental inclusion of unrelated demultiplexed data, and follows the input model used by established Nextflow pipelines such as `nf-core/atacseq`.

