# dmagna-atac

`dmagna-atac` is a Nextflow pipeline for processing paired-end ATAC-seq data from *Daphnia magna*. It takes raw FASTQ files, an explicit sample sheet, and a matching gene annotation as input. It produces filtered alignment files, CPM-normalized bigWig tracks, deepTools signal profiles and quality-control plots, per-sample peak calls, and run-level reports.

Developed as a KU Leuven research-course project, the pipeline runs on Nextflow and Apptainer. It takes inspiration from [nf-core/atacseq](https://github.com/nf-core/atacseq) and is an independent implementation with a focused initial scope.

## Current pipeline

The current implementation covers read processing through per-sample peak calling. Sample selection and metadata come from a CSV sample sheet, and the workflow supports a variable number of samples without changes to the code.

```mermaid
flowchart TD
    A["Sample sheet + paired-end FASTQ"] --> B["Validate input and group lanes by sample"]
    B --> C["Raw-read QC · FastQC"]
    B --> D["Merge lanes and trim adapters · Cutadapt"]
    D --> E["Trimmed-read QC · FastQC"]
    D --> F["Align reads · BWA-MEM2"]
    F --> G["Sort, deduplicate and filter BAM · Samtools"]
    G --> H["CPM bigWig coverage"]
    H --> I["deepTools · TSS and gene-body matrices, profiles and heatmaps"]
    G --> J["deepTools · plotFingerprint"]
    I --> K["Per-sample peak calling · MACS3 BAMPE"]
    C --> L["Aggregate QC reports · MultiQC"]
    D --> L
    E --> L
    G --> L
    I --> L
    J --> K
    K --> L
```

| Stage | Current implementation |
| --- | --- |
| Input validation | Checks required columns, sample IDs, file availability, gzip signatures, R1/R2 filename pairing, duplicate file entries, and metadata consistency. Saves the normalized sample sheet and validation report. |
| Lane handling | Combines read pairs sharing a `sample_id` before trimming, producing one analysis unit per sample. |
| Read QC and trimming | Runs FastQC before and after Cutadapt. Removes adapters and discards pairs if either read is shorter than the minimum length, which defaults to 20 bp. Additional Cutadapt options can be supplied with `--extra_cutadapt`. |
| Alignment | Aligns paired reads with BWA-MEM2 and produces sorted, indexed BAM files. |
| BAM processing | Removes duplicates, mitochondrial reads, unmapped or multi-mapped reads, pairs with more than four edit-distance mismatches per read, soft-clipped reads, fragments longer than 2 kb, interchromosomal pairs, and pairs outside FR orientation. Pair-level filters retain or discard both mates together. Blacklist filtering is not enabled. Produces indexed filtered BAMs, `flagstat` summaries, and a filter-count table. |
| Coverage and deepTools | Creates CPM-normalized bigWigs from filtered BAMs. `computeMatrix`, `plotProfile`, and `plotHeatmap` summarize signal around gene TSSs and across gene bodies. `plotFingerprint` uses the indexed filtered BAMs. The deepTools stage completes before peak calling. |
| Peak calling | Runs MACS3 independently for each sample in `BAMPE` mode using the supplied effective genome size. |
| Reporting | Aggregates QC reports with MultiQC and records software versions, validated inputs, filter metrics, and Nextflow execution reports. |

These are the stages implemented in the repository. Run completion and data quality are assessed from the corresponding VSC execution and QC reports.

## Inputs

The required inputs are:

| Parameter | Input |
| --- | --- |
| `--samplesheet` | CSV listing sample IDs and paired FASTQ paths. |
| `--reference_fasta` | Reference genome FASTA. |
| `--bwa_index_dir` | BWA-MEM2 index files generated from that FASTA, with the FASTA filename as their prefix. |
| `--gene_annotation` | GFF3 or GTF gene annotation whose sequence IDs and coordinates match the reference FASTA. |
| `--mitochondrial_contig` | Mitochondrial sequence name in the reference. |
| `--macs_genome_size` | Effective genome size used for peak calling. |
| `--container` | Path to the Apptainer SIF image. |

The sample sheet requires `sample_id`, `fastq_1`, and `fastq_2`. The optional `condition` and `replicate` columns carry experimental metadata; they do not currently trigger group comparisons or statistical testing.

```csv
sample_id,fastq_1,fastq_2,condition,replicate
sample_A,/path/to/sample_A_L001_R1.fastq.gz,/path/to/sample_A_L001_R2.fastq.gz,fish,1
sample_A,/path/to/sample_A_L002_R1.fastq.gz,/path/to/sample_A_L002_R2.fastq.gz,fish,1
sample_B,/path/to/sample_B_L001_R1.fastq.gz,/path/to/sample_B_L001_R2.fastq.gz,control,1
```

Use the same `sample_id` for multiple lanes of one sample and a different ID for each biological replicate. Repeated rows must have consistent metadata. Only files explicitly listed in the sample sheet are processed. Use absolute FASTQ paths accessible on the compute nodes.

The current project input selection contains `fish_r4`, `fish_r5`, `control_r4`, and `control_r1`; these names are not hard-coded into the workflow. A [sample sheet example](examples/samplesheet.example.csv) is included in the repository, while the actual sheet with VSC storage paths is maintained alongside the data.

## Outputs

Results are written under `--outdir` (default: `results/`):

| Directory | Main contents |
| --- | --- |
| `alignment/` | Sorted BAM files and indexes before filtering. |
| `filtered_bam/` | Final filtered BAM files, indexes, `flagstat` summaries, and pair-filter metrics. |
| `bigwig/` | CPM-normalized bigWig coverage tracks, with a default 10 bp bin size. |
| `deeptools/` | TSS and gene-body matrices, profiles, heatmaps, fingerprint plots, count tables, and the derived gene BED file. |
| `peaks/` | Per-sample MACS3 `narrowPeak`, summit BED, and peak-report files. |
| `qc/` | Published raw FastQC output, Cutadapt JSON reports, and the combined MultiQC report incorporating raw and trimmed-read QC. |
| `pipeline_info/` | Validated sample sheet, input-file manifest, validation report, software versions, and Nextflow execution report, timeline, trace, and DAG. |

The filtered BAMs, signal tracks, deepTools plots, and per-sample peak files provide the starting point for downstream analyses.

### Gene annotation for deepTools

The current project reference is `GCA_030254905.1` (`UOB_LRV0_1`). The matching assembly is [listed by NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030254905.1/), and a gene set for this isolate is described in [Chaturvedi et al. (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10570034/). Use a GFF3 or GTF from that exact assembly, and confirm its sequence IDs match the FASTA used to build the BWA-MEM2 index. The workflow accepts compressed `.gff.gz`, `.gff3.gz`, `.gtf.gz`, and uncompressed GFF3/GTF files. It creates gene-level intervals from `gene` features, or from transcript features if gene rows are absent. It does not substitute an annotation from another *D. magna* strain or assembly. Provide the selected file with `--gene_annotation`.

## Execution and reproducibility

- **Nextflow DSL2** manages dependencies, sample processing, and resumption with `-resume`.
- **Apptainer** supplies the software environment defined by `environment.yml` and `containers/dmagna-atac.def`.

See the [VSC execution guide](docs/vsc-run.md) for setup and submission details. The repository includes a [Slurm launcher](slurm/run-dmagna-atac.slurm); review its account, notification settings, module names, and storage paths for your environment before submission.

## Project reference configuration

The project uses the following *D. magna* reference settings, recorded in the existing project documentation:

| Field | Value |
| --- | --- |
| Assembly accession | `GCA_030254905.1` |
| Assembly name | `UOB_LRV0_1` |
| Genotype | `LRV0_1` |
| FASTA sequences | 72 |
| Total FASTA length | `129994773` bp |
| Mitochondrial contig | `CM058182.1` |
| MACS3 effective genome size | `129954900` |

The project defines the effective genome size as the number of non-`N` nuclear bases after excluding the mitochondrial sequence. Reference files, metadata, and generated indexes are maintained in the VSC data area.

## Planned additions

Future development may extend the current per-sample workflow toward comparisons between experimental conditions. The following features are planned and are not implemented in the current pipeline:

| Planned feature | Purpose |
| --- | --- |
| Additional ATAC-seq-specific QC | Add fragment-length distributions, TSS enrichment scores, and fraction of reads in peaks (FRiP). |
| Replicate assessment and consensus peaks | Evaluate agreement between replicates and construct a shared set of regions for comparisons across samples. |
| Peak quantification and differential accessibility | Build a fragment-count matrix over consensus regions and compare conditions using a confirmed replicate structure and statistical design. |
| Peak annotation | Associate accessible regions with genes and genomic features using annotations matched to the reference assembly. |

Development priorities will be refined after reviewing the VSC run results and sample quality. Group-level analyses will require confirmed experimental metadata and an appropriate replicate design.

## Repository structure

```text
main.nf                         Nextflow workflow
nextflow.config                 process resources and VSC execution profile
nextflow_schema.json            parameter schema
environment.yml                 pinned software environment
containers/dmagna-atac.def      Apptainer image definition
bin/validate_samplesheet.py     input validation and normalization
bin/filter_paired_alignments.py  pair-aware ATAC-seq alignment filtering
bin/make_deeptools_regions.py    gene-level BED generation for deepTools
examples/                       sample-sheet examples
docs/pipeline-design.md         workflow design and scope
docs/vsc-run.md                 VSC storage and execution guide
slurm/run-dmagna-atac.slurm      Slurm launcher for VSC runs
tests/                          validation tests
```

## Documentation

- [Pipeline design](docs/pipeline-design.md) — workflow contract, process boundaries, and Version 1 scope.
- [VSC execution guide](docs/vsc-run.md) — repository, storage, Apptainer, reference, and Slurm conventions.
- [Sample sheet example](examples/samplesheet.example.csv) — input format with placeholder paths.
