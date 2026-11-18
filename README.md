# dmagna_atac

An independent paired-end ATAC-seq pipeline for *Daphnia magna* implemented with Nextflow. Version 1 uses the sample sheet as the only source of input selection: it does not scan FASTQ directories or infer sample identity from filenames.

## Current implementation

- Strict validation and normalization of the CSV sample sheet before analysis starts.
- Multiple lanes can be assigned to the same `sample_id`; corresponding R1 and R2 files are merged before trimming.
- FastQC on raw FASTQ files.
- Adapter and quality trimming with Cutadapt.
- Paired-end alignment with BWA-MEM2, followed by BAM sorting and indexing.
- `samtools fixmate`, duplicate removal with `markdup`, proper-pair/alignment filtering, and mitochondrial-contig removal.
- Per-sample peak calling with MACS3 in `BAMPE` mode.
- QC for trimmed reads, MultiQC, and Nextflow execution reports, timeline, and trace.
- Recording of `samplesheet.used.csv`, `input_files.tsv`, the input-validation report, and software versions.

The reference FASTA, BWA-MEM2 index, mitochondrial contig name, and MACS3 effective genome size must be explicitly provided for the selected *D. magna* reference. They are intentionally not hard-coded in the pipeline.

## Run the sample-sheet tests

```bash
python3 -m unittest discover -s tests -v
```

## Build the container

The repository uses an Apptainer definition file as the container source:

```text
containers/dmagna-atac.def
```

Build the versioned SIF image on VSC, preferably from a compute allocation with enough memory:

```bash
apptainer build --fakeroot \
  "$VSC_DATA/dmagna_atac/containers/dmagna-atac-0.1.0.sif" \
  "$HOME/dmagna_atac/containers/dmagna-atac.def"
```

The definition installs the exact versions listed in `environment.yml` and runs built-in container tests during the build. The SIF image is generated on VSC and is not stored in Git.

## Prepare the reference genome

Select and record the final *D. magna* reference version, then prepare the FASTA and BWA-MEM2 index. The following is an example:

```bash
bwa-mem2 index -p /path/to/reference/bwa_index/dmagna.fa /path/to/reference/dmagna.fa
```

`--bwa_index_dir` must point to the directory containing all index files generated from that prefix. The mitochondrial contig name and `--macs_genome_size` must come from the same reference configuration; neither should be inferred from sample filenames.

## Run on VSC

The complete VSC procedure is documented in [`docs/vsc-run.md`](docs/vsc-run.md). The short form is:

```bash
nextflow run . \
  -profile vsc \
  -resume \
  -work-dir "$VSC_SCRATCH/dmagna_atac/work" \
  --samplesheet "$VSC_DATA/dmagna_atac/metadata/samplesheet.csv" \
  --reference_fasta "$VSC_DATA/dmagna_atac/reference/dmagna.fa" \
  --bwa_index_dir "$VSC_DATA/dmagna_atac/reference/bwa_index" \
  --mitochondrial_contig '<MITO_CONTIG_FROM_REFERENCE>' \
  --macs_genome_size '<EFFECTIVE_GENOME_SIZE>' \
  --container "$VSC_DATA/dmagna_atac/containers/dmagna-atac-0.1.0.sif" \
  --outdir "$VSC_DATA/dmagna_atac/results/run-001"
```

The `vsc` profile submits processes to Slurm and runs them with the local Apptainer SIF image. If the VSC allocation requires explicit values, add `--slurm_account <account>` and `--slurm_partition <partition>`.

Before the first full run, use a small real paired-end subset to verify the container, reference index, mitochondrial contig, and scheduler configuration. Preserve the input manifest and execution records under `results/pipeline_info/` for every formal run.

## Output structure

The main outputs are organized as follows:

```text
results/
├── qc/raw/<sample_id>/
├── qc/trimmed/<sample_id>/
├── qc/multiqc/
├── alignment/<sample_id>.sorted.bam[.bai]
├── filtered_bam/<sample_id>.filtered.bam[.bai]
├── peaks/<sample_id>.narrowPeak
└── pipeline_info/
```

## Next steps

Before the formal analysis, complete these tasks:

1. Pull the repository on VSC and build `dmagna-atac-0.1.0.sif` from the Apptainer definition.
2. Confirm `condition` and `replicate` for the four samples from the experimental records, then add the real FASTQ paths to the VSC sample sheet.
3. Confirm the reference version, mitochondrial contig name, and MACS3 effective genome size; generate and verify the BWA-MEM2 index on VSC.
4. Run a small VSC smoke test, then run the complete four-sample workflow. Inspect `input_validation_report.txt`, MultiQC, BAM statistics, and the peak files for every sample.

Differential accessibility, consensus peaks, peak annotation, motif analysis, and functional enrichment are outside the scope of Version 1. They should be designed separately after the four-sample QC and biological-replicate structure have been confirmed.
