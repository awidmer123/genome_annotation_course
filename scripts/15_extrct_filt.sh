#!/bin/bash
#SBATCH --job-name=Extract_and_filtering
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/15_%x_%j.out
#SBATCH --error=logs/15_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

# setting paths
WORKDIR=/data/users/awidmer/genome_annotation_course
OUTDIR="$WORKDIR/output/09_MAKER"
FINALDIR="$OUTDIR/final"
COURSEDIR=/data/courses/genome_annotation_resources/2026
MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
ANALYSIS_SIF=$COURSEDIR/containers/Annotation_analysis_R4.4.sif

set -euo pipefail

gff=assembly.all.maker.noseq.gff
protein=assembly.all.maker.proteins.fasta
transcript=assembly.all.maker.transcripts.fasta

mkdir -p "$OUTDIR/tmp"
cd "$FINALDIR"

MAKER="apptainer exec --cleanenv --env TMPDIR=/tmp --bind $OUTDIR/tmp:/tmp --bind $COURSEDIR --bind $WORKDIR --pwd $FINALDIR $MAKER_SIF"
ANALYSIS="apptainer exec --cleanenv --bind $COURSEDIR --bind $WORKDIR --pwd $FINALDIR $ANALYSIS_SIF"

# 4. AED-Verteilung
$MAKER perl /usr/local/bin/AED_cdf_generator.pl -b 0.025 ${gff}.renamed.gff \
    > assembly.all.maker.renamed.gff.AED.txt

# Bonus für die Guiding Question: Pfam-Anteil bei tiefem vs. hohem AED
awk -F'\t' '$3=="mRNA" {
    match($9, /_AED=[0-9.]+/); aed = substr($9, RSTART+5, RLENGTH-5) + 0
    bin = (aed <= 0.5) ? "AED<=0.5" : "AED>0.5"
    n[bin]++; if ($9 ~ /Pfam:/) p[bin]++; tot++; ptot += ($9 ~ /Pfam:/)
} END {
    for (b in n) printf "%s: %d mRNAs, %d mit Pfam (%.1f%%)\n", b, n[b], p[b], 100*p[b]/n[b]
    printf "Total: %d mRNAs, %d mit Pfam (%.1f%%)\n", tot, ptot, 100*ptot/tot
}' ${gff}.renamed.iprscan.gff > pfam_by_AED.txt

# 5. Quality filter (AED < 1 und/oder Pfam-Domäne)
$MAKER perl /usr/local/bin/quality_filter.pl -s ${gff}.renamed.iprscan.gff \
    > ${gff}_iprscan_quality_filtered.gff

# 6. Nur Gen-Features behalten
grep -P "\tgene\t|\tCDS\t|\texon\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t|\tmRNA\t" \
    ${gff}_iprscan_quality_filtered.gff > filtered.genes.renamed.gff3
sed -i '1i##gff-version 3' filtered.genes.renamed.gff3
echo "Feature types:"; cut -f3 filtered.genes.renamed.gff3 | grep -v '^#' | sort | uniq -c

# 7. Gefilterte Transkripte/Proteine extrahieren
grep -P "\tmRNA\t" filtered.genes.renamed.gff3 | awk '{print $9}' | cut -d ';' -f1 | sed 's/ID=//g' > list.txt
$ANALYSIS seqkit grep -f list.txt ${transcript}.renamed.fasta -o ${transcript}.renamed.filtered.fasta
$ANALYSIS seqkit grep -f list.txt ${protein}.renamed.fasta    -o ${protein}.renamed.filtered.fasta

echo "mRNAs nach Filter: $(wc -l < list.txt)"
echo "Proteine: $(grep -c '>' ${protein}.renamed.filtered.fasta)"