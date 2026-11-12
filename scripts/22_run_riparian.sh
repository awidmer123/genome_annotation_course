#!/usr/bin/env bash
#SBATCH --job-name=riparian
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/22_%x_%j.out
#SBATCH --error=logs/22_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

set -euo pipefail

# setting paths
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
FINALDIR="$WORKDIR/output/09_MAKER/final"
ACCESSION=Mh_0
OTHERS=(Taz_0 Altai_5 Ice_1)

PROTEINS="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"
LONGEST_PROTEINS="$FINALDIR/assembly.all.maker.proteins.fasta.renamed.longest.fasta"
GENES="$FINALDIR/filtered.genes.renamed.gff3"

UNIPROT="$COURSEDIR/data/uniprot_viridiplantae_reviewed.fa"
TAIR10="$COURSEDIR/data/TAIR10_pep_20110103_representative_gene_model"

OUTDIR="$WORKDIR/output/comparisons"
GENESPACE_DIR="$OUTDIR/GENESPACE"

MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
GENESPACE_SIF="$COURSEDIR/containers/genespace_latest.sif"
ANALYSIS_SIF="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"
# -----------------------------------------

# Scaffolds that seemed flipped
FLIPS=()

TMPDIR_HOST="$OUTDIR/GENESPACE_R_tmp"
mkdir -p "$TMPDIR_HOST"
apptainer exec --cleanenv \
    --bind /data --bind "$WORKDIR" \
    --bind "$TMPDIR_HOST":/temp --env TMPDIR=/temp --pwd "$WORKDIR" \
    "$GENESPACE_SIF" Rscript "$WORKDIR/scripts/22_riparian.R" "$GENESPACE_DIR" "$ACCESSION" "${FLIPS[@]}"
