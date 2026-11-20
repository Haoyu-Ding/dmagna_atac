# Running `dmagna_atac` on KU Leuven VSC

This document describes the first formal execution on VSC. The local computer is used for editing and Git synchronization only. FASTQ processing is submitted to the VSC Slurm scheduler.

## 1. Repository and storage layout

Keep the repository in the home directory and keep large data outside Git:

```text
$HOME/dmagna_atac/
    workflow code and configuration

$VSC_DATA/dmagna_atac/
    reference/
    metadata/
    containers/
    results/

$VSC_SCRATCH/dmagna_atac/
    reads/
    work/
    logs/
```

Create the directories once:

```bash
mkdir -p "$VSC_DATA/dmagna_atac"/{reference,metadata,containers,results}
mkdir -p "$VSC_SCRATCH/dmagna_atac"/{reads,work,logs,apptainer-tmp}
```

The `reference`, `metadata`, and `results` directories are long-term data. The `reads` and `work` directories are scratch data and must not be treated as backups.

## 2. Clone or update the repository

The repository is cloned into `$HOME`, not into the data or scratch directory:

```bash
cd "$HOME"
git clone https://github.com/Haoyu-Ding/dmagna_atac.git
cd dmagna_atac
```

For later updates:

```bash
cd "$HOME/dmagna_atac"
git pull --ff-only
```

## 3. Load the VSC software modules

Module names can differ between VSC clusters and toolchains. Check the available names first, then load one Nextflow module and one Apptainer module:

```bash
module spider Nextflow
module spider Apptainer
module load Nextflow/<available-version>
module load Apptainer/<available-version>
nextflow -version
apptainer --version
```

Do not put a guessed module version in the repository. Record the actual modules used in the run notes.

## 4. Prepare the Apptainer environment

The cache is kept in `VSC_DATA`, while temporary image conversion files are kept in `VSC_SCRATCH`:

```bash
export APPTAINER_CACHEDIR="$VSC_DATA/dmagna_atac/containers/cache"
export APPTAINER_TMPDIR="$VSC_SCRATCH/dmagna_atac/apptainer-tmp"
export NXF_APPTAINER_CACHEDIR="$VSC_DATA/dmagna_atac/containers/cache"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
```

The container source is the Apptainer definition committed in the repository:

```text
$HOME/dmagna_atac/containers/dmagna-atac.def
```

Build the versioned SIF image on a compute node. The exact account or partition flags depend on the VSC allocation; append them if required by your project:

```bash
cd "$HOME/dmagna_atac"

srun --time=01:00:00 --mem=16G --cpus-per-task=4 \
  apptainer build --fakeroot \
  "$VSC_DATA/dmagna_atac/containers/dmagna-atac-0.1.0.sif" \
  "$HOME/dmagna_atac/containers/dmagna-atac.def"
```

The `%test` section of the definition checks Python, FastQC, Cutadapt, BWA-MEM2, Samtools, MACS3, MultiQC, and the sample-sheet validator during the build. If the cluster does not allow `--fakeroot`, use the VSC-supported Apptainer build method for your account or contact VSC support; do not replace the definition with a separate unrecorded software installation.

## 5. Transfer FASTQ files with Globus

Transfer only the selected paired-end reads to:

```text
$VSC_SCRATCH/dmagna_atac/reads/
```

For the current run, include only the four selected samples and their `R1`/`R2` files. Do not include `I1`, `I2`, `Undetermined`, or the two excluded barcode libraries.

Keep the FASTQ files compressed as `.fastq.gz`. Do not add them to Git.

## 6. Create the actual sample sheet on VSC

Create `$VSC_DATA/dmagna_atac/metadata/samplesheet.csv` using the real VSC paths:

```csv
sample_id,fastq_1,fastq_2,condition,replicate
fish_r4,/scratch/.../GC187803_CCTGTCAT-ATGACAGG_S17_L001_R1_001.fastq.gz,/scratch/.../GC187803_CCTGTCAT-ATGACAGG_S17_L001_R2_001.fastq.gz,fish,
fish_r5,/scratch/.../GC187803_TTAGGTCG-CGACCTAA_S22_L001_R1_001.fastq.gz,/scratch/.../GC187803_TTAGGTCG-CGACCTAA_S22_L001_R2_001.fastq.gz,fish,
control_r4,/scratch/.../GC187803_GTGCCATA-TATGGCAC_S20_L001_R1_001.fastq.gz,/scratch/.../GC187803_GTGCCATA-TATGGCAC_S20_L001_R2_001.fastq.gz,control,
control_r1,/scratch/.../GC187803_GGCGTTAT-ATAACGCC_S19_L001_R1_001.fastq.gz,/scratch/.../GC187803_GGCGTTAT-ATAACGCC_S19_L001_R2_001.fastq.gz,control,
```

The `condition` and `replicate` values must be confirmed from the experimental record. They are metadata only in Version 1.

## 7. Prepare the reference

Place the selected reference FASTA and its BWA-MEM2 index in the same reference area:

```text
$VSC_DATA/dmagna_atac/reference/GCA_030254905.1/dmagna_GCA_030254905.1.fa
$VSC_DATA/dmagna_atac/reference/GCA_030254905.1/index/dmagna_GCA_030254905.1.fa.*
```

The selected assembly is `GCA_030254905.1` (`UOB_LRV0_1`, genotype `LRV0_1`). The index files must be generated from this exact FASTA. For this assembly, use `CM058182.1` as the mitochondrial contig and `129954900` as the effective genome size. These values are required command-line parameters and are deliberately not guessed by the workflow.

## 8. Submit the pipeline

Launch Nextflow from a VSC login node, where `sbatch` is available. Nextflow submits the individual processes to Slurm. Use a persistent terminal session if the SSH connection may close.

Example with a local SIF image:

```bash
cd "$HOME/dmagna_atac"

nextflow run . \
  -profile vsc \
  -resume \
  -work-dir "$VSC_SCRATCH/dmagna_atac/work" \
  --samplesheet "$VSC_DATA/dmagna_atac/metadata/samplesheet.csv" \
  --reference_fasta "$VSC_DATA/dmagna_atac/reference/GCA_030254905.1/dmagna_GCA_030254905.1.fa" \
  --bwa_index_dir "$VSC_DATA/dmagna_atac/reference/GCA_030254905.1/index" \
  --mitochondrial_contig 'CM058182.1' \
  --macs_genome_size '129954900' \
  --container "$VSC_DATA/dmagna_atac/containers/dmagna-atac-0.1.0.sif" \
  --outdir "$VSC_DATA/dmagna_atac/results/run-001"
```

If the project allocation requires an account or a non-default partition, append:

```text
--slurm_account <project-account> --slurm_partition <partition-name>
```

Do not fill these values from memory; use the allocation information provided for your VSC project.

## 9. Verify a completed run

The run is complete only when Nextflow reports success and the result directory contains:

```text
results/qc/multiqc/multiqc_report.html
results/filtered_bam/<sample_id>.filtered.bam
results/peaks/<sample_id>.narrowPeak
results/pipeline_info/samplesheet.used.csv
results/pipeline_info/input_validation_report.txt
results/pipeline_info/software_versions.yml
```

Inspect the validation report, MultiQC report, BAM flagstat files, and peak files for all samples before treating the run as biologically interpretable.
