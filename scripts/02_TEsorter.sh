#!/bin/bash
#SBATCH --job-name=TEsorter
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --reservation=class-sbl-30004-2026-10-14
#SBATCH --output=logs/02_%x_%j.out
#SBATCH --error=logs/02_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR=/data/users/awidmer/genome_annotation_course  # Directory in your /data/users/$USER/Folder.
# Use the intact-element FASTA
# Edit these two input paths to match your assembly's filename.
RAW="$WORKDIR/output/01_EDTA"
# Expand the wildcard outside quotes;
shopt -s nullglob
FASTA=("$RAW"/*.mod.EDTA.raw/*.mod.LTR.intact.raw.fa)
GFF=("$RAW"/*.mod.EDTA.raw/*.mod.LTR.intact.raw.gff3)

OUTDIR="$WORKDIR/output/02_LTR_full_length"
unset SINGULARITY_BIND APPTAINER_BIND SINGULARITY_TMPDIR
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/tmp"

# Work with copies so ongoing EDTA files are not changed.
cp "${FASTA[0]}" "$OUTDIR/intact_LTR.fa"
cp "${GFF[0]}" "$OUTDIR/intact_LTR.gff3"
cd "$OUTDIR"


apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR" \
    --bind "$OUTDIR/tmp":/tmp --pwd "$OUTDIR" \
    "$COURSEDIR/containers/TEsorter_1.3.0.sif" \
    TEsorter intact_LTR.fa -db rexdb-plant -p "$SLURM_CPUS_PER_TASK"

# After this job succeeds, run the R script (04_annotation_circlize.R) in the output directory
mkdir -p plots
cd "$WORKDIR/output/02_LTR_full_length"
apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR" --pwd "$WORKDIR/output/02_LTR_full_length" \
"$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
Rscript --vanilla "$WORKDIR/scripts/03_LTR_identity.R"