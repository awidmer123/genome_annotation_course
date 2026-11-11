#!/usr/bin/env bash
#SBATCH --job-name=gs_prep
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/18_%x_%j.out
#SBATCH --error=logs/18_%x_%j.err
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

mkdir -p "$GENESPACE_DIR/bed" "$GENESPACE_DIR/peptide" "$GENESPACE_DIR/tmp"

# 1. Eigene Accession: BED (0-based) aus Gene-Features
grep -P "\tgene\t" "$GENES" \
  | awk 'BEGIN{OFS="\t"} {split($9,a,";"); split(a[1],b,"="); print $1, $4-1, $5, b[2]}' \
  > "$GENESPACE_DIR/bed/${ACCESSION}.bed"

# Longest-Protein-FASTA, Header auf Gen-ID kürzen (-RA usw. weg)
sed '/^>/s/ .*//' "$LONGEST_PROTEINS" | sed '/^>/s/-R[A-Z]*$//' \
  > "$GENESPACE_DIR/peptide/${ACCESSION}.fa"

# 2. TAIR10
cp "$COURSEDIR/data/TAIR10.bed" "$GENESPACE_DIR/bed/TAIR10.bed"
cp "$COURSEDIR/data/TAIR10.fa"  "$GENESPACE_DIR/peptide/TAIR10.fa"

# 3. Drei andere Accessions
for other in "${OTHERS[@]}"; do
    [[ "$other" == "$ACCESSION" ]] && { echo "FEHLER: $other ist deine eigene Accession"; exit 1; }
    cp "$COURSEDIR/data/Lian_et_al/GENESPACE/bed/${other}.bed"     "$GENESPACE_DIR/bed/"
    cp "$COURSEDIR/data/Lian_et_al/GENESPACE/peptide/${other}.fa"  "$GENESPACE_DIR/peptide/"
done

# 4. Validierung: Duplikate und BED <-> FASTA Abgleich
cd "$GENESPACE_DIR"
fail=0
for bed in bed/*.bed; do
    g=$(basename "$bed" .bed)
    fa="peptide/${g}.fa"
    cut -f4 "$bed" | sort > tmp/bed_ids
    grep '>' "$fa" | sed 's/^>//' | sort > tmp/fa_ids
    dup_bed=$(uniq -d tmp/bed_ids | wc -l)
    dup_fa=$(uniq -d tmp/fa_ids | wc -l)
    missing=$(comm -23 <(sort -u tmp/bed_ids) <(sort -u tmp/fa_ids) | wc -l)   # in BED, nicht in FASTA
    extra=$(comm -13 <(sort -u tmp/bed_ids) <(sort -u tmp/fa_ids) | wc -l)     # in FASTA, nicht in BED
    printf "%-12s genes=%6d  dupBED=%d dupFA=%d missing=%d extra=%d\n" \
        "$g" "$(wc -l < tmp/bed_ids)" "$dup_bed" "$dup_fa" "$missing" "$extra"
    (( dup_bed + dup_fa + missing + extra > 0 )) && fail=1
done
(( fail )) && { echo "ID-Probleme gefunden, erst fixen!"; exit 1; }
echo "Alle 5 Paare OK."
