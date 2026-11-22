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

This project targets the KU Leuven/UHasselt `wICE` cluster. Check the available module names first, then load explicit versions:

```bash
module spider Nextflow
module spider Apptainer
module load Nextflow/<available-version>
module load Apptainer/<available-version>
nextflow -version
apptainer --version
```

Do not put a guessed module version in the repository. Record the actual modules used in the run notes. VSC job scripts should use a login shell (`#!/bin/bash -l`) and load the required modules inside the job environment.

## 4. Prepare the Apptainer environment

For a one-off container build, use node-local temporary directories as recommended by VSC. The immutable SIF image itself is stored in `VSC_DATA`:

```bash
export APPTAINER_CACHEDIR="/tmp/$USER/apptainer_cachedir"
export APPTAINER_TMPDIR="/tmp/$USER/apptainer_tmpdir"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
```

The container source is the Apptainer definition committed in the repository:

```text
$HOME/dmagna_atac/containers/dmagna-atac.def
```

Check which Slurm credit account is available to you:

```bash
sam-balance
```

Choose the account shown by that command. Every KU Leuven Tier-2 job must specify a valid account. The `interactive` partition is appropriate for this one-off container build and allows up to eight cores.

Start an interactive wICE compute allocation, replacing `YOUR_ACCOUNT` with the account returned by `sam-balance`:

```bash
ACCOUNT='YOUR_ACCOUNT'

srun --account="$ACCOUNT" \
  --clusters=wice \
  --partition=interactive \
  --ntasks=1 \
  --cpus-per-task=4 \
  --mem=16G \
  --time=01:00:00 \
  --pty bash -l
```

Run the build inside the allocated shell:

```bash
cd "$HOME/dmagna_atac"
module load Apptainer/<available-version>

export APPTAINER_CACHEDIR="/tmp/$USER/apptainer_cachedir"
export APPTAINER_TMPDIR="/tmp/$USER/apptainer_tmpdir"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

apptainer build --fakeroot \
  "$VSC_DATA/dmagna_atac/containers/dmagna-atac-0.1.1.sif" \
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

Submit the project Slurm launcher from a VSC login node. The launcher itself runs Nextflow in a small `batch` job, while the `vsc` profile submits every workflow process to Slurm. This keeps the run independent of the SSH connection and uses the local Apptainer SIF image.

```bash
cd "$HOME/dmagna_atac"

module spider Nextflow
module spider Apptainer

export NEXTFLOW_MODULE='Nextflow/<available-version>'
export APPTAINER_MODULE='Apptainer/<available-version>'

# Edit the account line locally before submitting. Use the value returned by
# sam-balance, and do not commit your real account to the public repository.
${EDITOR:-vi} "$HOME/dmagna_atac/slurm/run-dmagna-atac.slurm"

mkdir -p "$VSC_DATA/dmagna_atac/logs"

sbatch \
  --clusters=wice \
  --output="$VSC_DATA/dmagna_atac/logs/dmagna-atac-%j.out" \
  --error="$VSC_DATA/dmagna_atac/logs/dmagna-atac-%j.err" \
  --export=ALL,NEXTFLOW_MODULE,APPTAINER_MODULE \
  "$HOME/dmagna_atac/slurm/run-dmagna-atac.slurm"
```

Replace both module placeholders with names shown by `module spider`. In the Slurm file, replace `#SBATCH --account=YOUR_ACCOUNT` with a valid account from `sam-balance`. The Slurm file supplies the reference, samplesheet, SIF, work, and output paths defined for this project. To run a new analysis, override `OUTDIR` and, when appropriate, `WORK_DIR` in the submission environment.

Monitor the launcher and its child tasks with:

```bash
squeue --clusters=wice --user "$USER"
tail -f "$VSC_DATA/dmagna_atac/logs/dmagna-atac-<JOBID>.out"
```

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
