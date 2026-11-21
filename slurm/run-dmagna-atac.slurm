#!/bin/bash -l
#SBATCH --clusters=wice
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=batch
#SBATCH --job-name=dmagna-atac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

set -euo pipefail

die() {
    echo "ERROR: $*" >&2
    exit 1
}

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/dmagna_atac}"
DATA_ROOT="${VSC_DATA:?VSC_DATA is not set}"
SCRATCH_ROOT="${VSC_SCRATCH:?VSC_SCRATCH is not set}"
ACCOUNT="${SLURM_JOB_ACCOUNT:-}"
PIPELINE_PARTITION="${SLURM_JOB_PARTITION:-batch}"

[[ -n "$ACCOUNT" ]] || die "No Slurm account found; submit with --account=<account>"
[[ "$ACCOUNT" != "YOUR_ACCOUNT" ]] || die "Replace YOUR_ACCOUNT in the #SBATCH --account line before submission"

SIF="${SIF:-$DATA_ROOT/dmagna_atac/containers/dmagna-atac-0.1.0.sif}"
SAMPLESHEET="${SAMPLESHEET:-$DATA_ROOT/dmagna_atac/metadata/samplesheet.csv}"
REFERENCE_FASTA="${REFERENCE_FASTA:-$DATA_ROOT/dmagna_atac/reference/GCA_030254905.1/dmagna_GCA_030254905.1.fa}"
BWA_INDEX_DIR="${BWA_INDEX_DIR:-$DATA_ROOT/dmagna_atac/reference/GCA_030254905.1/index}"
OUTDIR="${OUTDIR:-$DATA_ROOT/dmagna_atac/results/run-001}"
WORK_DIR="${WORK_DIR:-$SCRATCH_ROOT/dmagna_atac/work}"

if ! command -v nextflow >/dev/null 2>&1; then
    : "${NEXTFLOW_MODULE:?Nextflow is not available; set NEXTFLOW_MODULE to the VSC module name}"
    module load "$NEXTFLOW_MODULE"
fi

if ! command -v apptainer >/dev/null 2>&1; then
    : "${APPTAINER_MODULE:?Apptainer is not available; set APPTAINER_MODULE to the VSC module name}"
    module load "$APPTAINER_MODULE"
fi

[[ -d "$PROJECT_ROOT" ]] || die "Project directory not found: $PROJECT_ROOT"
[[ -r "$SIF" ]] || die "Container not readable: $SIF"
[[ -r "$SAMPLESHEET" ]] || die "Samplesheet not readable: $SAMPLESHEET"
[[ -r "$REFERENCE_FASTA" ]] || die "Reference FASTA not readable: $REFERENCE_FASTA"
[[ -d "$BWA_INDEX_DIR" ]] || die "BWA index directory not found: $BWA_INDEX_DIR"

INDEX_PREFIX="$BWA_INDEX_DIR/$(basename "$REFERENCE_FASTA")"
compgen -G "${INDEX_PREFIX}.*" >/dev/null || die "No BWA index files found for $INDEX_PREFIX"

mkdir -p "$OUTDIR" "$WORK_DIR"

export NXF_HOME="${NXF_HOME:-$DATA_ROOT/dmagna_atac/nextflow}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$SCRATCH_ROOT/dmagna_atac/apptainer-cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$SCRATCH_ROOT/dmagna_atac/apptainer-tmp}"
mkdir -p "$NXF_HOME" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "Starting dmagna_atac"
echo "Project: $PROJECT_ROOT"
echo "Container: $SIF"
echo "Account: $ACCOUNT"
echo "Pipeline partition: $PIPELINE_PARTITION"
echo "Output: $OUTDIR"

exec nextflow run "$PROJECT_ROOT" \
    -profile vsc \
    -resume \
    -work-dir "$WORK_DIR" \
    --samplesheet "$SAMPLESHEET" \
    --reference_fasta "$REFERENCE_FASTA" \
    --bwa_index_dir "$BWA_INDEX_DIR" \
    --mitochondrial_contig 'CM058182.1' \
    --macs_genome_size '129954900' \
    --container "$SIF" \
    --outdir "$OUTDIR" \
    --slurm_account "$ACCOUNT" \
    --slurm_partition "$PIPELINE_PARTITION"
