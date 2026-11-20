# dmagna-atac

An independent paired-end ATAC-seq workflow for *Daphnia magna*, developed as a KU Leuven research-course project.

The project takes inspiration from [nf-core/atacseq](https://github.com/nf-core/atacseq), but is implemented as a smaller, independently designed workflow. Its purpose is to provide a clear and reproducible starting point for learning workflow development, containerized software management, and ATAC-seq data processing on the KU Leuven VSC.

## Project goal

Version 1 focuses on a complete per-sample path from raw paired-end reads to reproducible peak calls. The workflow is designed to accept an arbitrary number of samples through an explicit sample sheet; sample identity is never inferred from directory contents or sequencing filenames.

The current project dataset contains four selected libraries:

| Sample | Condition | Source library |
| --- | --- | --- |
| `fish_r4` | fish | `S17`, `CCTGTCAT-ATGACAGG` |
| `fish_r5` | fish | `S22`, `TTAGGTCG-CGACCTAA` |
| `control_r4` | control | `S20`, `GTGCCATA-TATGGCAC` |
| `control_r1` | control | `S19`, `GGCGTTAT-ATAACGCC` |

Index reads, undetermined reads, and the two excluded libraries are not part of this Version 1 input set. The final sample sheet remains a VSC-side metadata file because it contains storage-specific FASTQ paths.

## Version 1 workflow

```text
paired-end FASTQ
        |
sample-sheet validation
        |
raw-read QC (FastQC)
        |
adapter and quality trimming (Cutadapt)
        |
paired-end alignment (BWA-MEM2)
        |
BAM sorting, indexing, duplicate removal, filtering
        |
mitochondrial-read removal
        |
per-sample peak calling (MACS3, BAMPE)
        |
QC aggregation and execution reports (MultiQC / Nextflow)
```

The main per-sample outputs are:

- filtered, indexed BAM files;
- MACS3 `narrowPeak` and summit files;
- raw and trimmed-read QC reports;
- alignment and filtering statistics;
- a combined MultiQC report;
- input, software, and workflow execution records.

Differential accessibility, consensus peak construction, peak annotation, motif analysis, and functional enrichment are intentionally outside the scope of Version 1.

## Reproducibility design

- **Nextflow** controls process dependencies, data flow, parallel execution, and resumability.
- **Apptainer** provides the versioned software environment required on VSC.
- **Slurm** is used for scheduled execution on VSC.
- **The sample sheet is the input contract.** It defines sample IDs, paired FASTQ files, and optional metadata such as condition and replicate.
- **Reference parameters are explicit.** The workflow requires the reference FASTA, BWA-MEM2 index, mitochondrial contig, and MACS3 effective genome size to be supplied together.
- **Large data are kept outside Git.** Code and configuration belong in the repository; FASTQ files, reference files, containers, intermediate files, and results belong in VSC storage.

## Current reference configuration

The selected working reference is the NCBI *D. magna* GenBank assembly:

| Field | Value |
| --- | --- |
| Assembly accession | `GCA_030254905.1` |
| Assembly name | `UOB_LRV0_1` |
| Genotype | `LRV0_1` |
| FASTA sequences | 72 |
| Total FASTA length | `129994773` bp |
| Mitochondrial contig | `CM058182.1` |
| MACS3 effective genome size | `129954900` |

The effective genome size is defined as non-`N` nuclear bases after excluding the mitochondrial sequence. The reference metadata and generated index are maintained in the VSC data area rather than committed to Git. The official assembly record is available from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030254905.1/).

## Repository structure

```text
main.nf                         Nextflow workflow
nextflow.config                 local and VSC execution profiles
nextflow_schema.json            parameter schema
environment.yml                 pinned software environment
containers/dmagna-atac.def      Apptainer image definition
bin/validate_samplesheet.py     input validation and normalization
examples/                       sample-sheet examples
docs/pipeline-design.md         workflow design and scope
docs/vsc-run.md                 VSC storage and execution guide
tests/                          validation tests
```

## Documentation

- [Pipeline design](docs/pipeline-design.md) — workflow contract, process boundaries, and Version 1 scope.
- [VSC execution guide](docs/vsc-run.md) — repository, storage, Apptainer, reference, and Slurm conventions.
- [Sample sheet example](examples/samplesheet.example.csv) — input format with placeholder paths.

## Project status

The workflow structure, sample-sheet validation, Apptainer definition, selected input libraries, and reference configuration are prepared. The remaining work is the VSC compute-stage verification: building the Apptainer image, generating the BWA-MEM2 index, running a smoke test, and evaluating the first formal four-sample run.
