#!/usr/bin/env bash
#SBATCH --job-name=pangenome
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/21_%x_%j.out
#SBATCH --error=logs/21_%x_%j.err
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

# Kursskript holen (einmalig) und wd / focal_genome eintragen
LOCAL_R="$WORKDIR/scripts/21-process_pangenome.R"
if [[ ! -f "$LOCAL_R" ]]; then
    SRC=$(ls "$COURSEDIR"/scripts/*process_pangenome.R | head -n1)
    cp "$SRC" "$LOCAL_R"
fi
sed -i -E "s|^(\s*wd\s*<-).*|\1 \"$GENESPACE_DIR\"|; s|^(\s*focal_genome\s*<-).*|\1 \"$ACCESSION\"|" "$LOCAL_R"
echo "Gesetzt:"; grep -nE '^\s*(wd|focal_genome)\s*<-' "$LOCAL_R"

TMPDIR_HOST="$OUTDIR/Annotation_R_tmp"
mkdir -p "$TMPDIR_HOST"
apptainer exec --cleanenv \
    --bind /data --bind "$WORKDIR" \
    --bind "$TMPDIR_HOST":/temp --env TMPDIR=/temp --pwd "$WORKDIR" \
    "$ANALYSIS_SIF" Rscript "$LOCAL_R"
