#!/usr/bin/env bash
#SBATCH --job-name=blastp
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/17_%x_%j.out
#SBATCH --error=logs/17_%x_%j.err
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

mkdir -p "$OUTDIR/tmp"
cd "$OUTDIR"
module load BLAST+/2.15.0-gompi-2021a

# 1. Both searches in parallel (DBs are already indexed); skipped if results exist
if [[ -s blastp_uniprot.tsv && -s blastp_TAIR10.tsv ]]; then
    echo "BLAST results found, skipping searches"
else
    blastp -query "$PROTEINS" -db "$UNIPROT" -num_threads 10 -outfmt 6 \
        -evalue 1e-5 -max_target_seqs 10 -out blastp_uniprot.tsv &
    uniprot_job=$!
    blastp -query "$PROTEINS" -db "$TAIR10" -num_threads 10 -outfmt 6 \
        -evalue 1e-5 -max_target_seqs 10 -out blastp_TAIR10.tsv &
    tair_job=$!
    wait "$uniprot_job"
    wait "$tair_job"
fi

# 2. Best hit per query
for ref in uniprot TAIR10; do
    sort -k1,1 -k12,12gr -k11,11g "blastp_${ref}.tsv" | awk -F '\t' '!seen[$1]++' \
        > "blastp_${ref}.besthits.tsv"
done

# 2b. BLAST+ 2.15 reports bare UniProt accessions (Q9LX89), but maker_functional_*
#     look up the full ID (sp|Q9LX89|FB255_ARATH). Without remapping,
#     maker_functional_fasta silently writes "Protein of unknown function" everywhere
#     and maker_functional_gff crashes ("Can't use string ("") as a HASH ref").
awk -F'\t' -v OFS='\t' '
  NR==FNR { if (/^>/) { id=substr($1,2); split(id,a,"|"); full[a[2]]=id }; next }
  { if ($2 in full) $2=full[$2]; print }' \
  <(cut -d' ' -f1 "$UNIPROT") blastp_uniprot.besthits.tsv \
  > blastp_uniprot.besthits.fullid.tsv

# 3. UniProt functions to gff3 and fasta
MAKER="apptainer exec --cleanenv --env TMPDIR=/tmp --bind /data --bind $WORKDIR --bind $OUTDIR/tmp:/tmp --pwd $OUTDIR $MAKER_SIF"
$MAKER maker_functional_fasta "$UNIPROT" blastp_uniprot.besthits.fullid.tsv "$PROTEINS" \
    > maker_proteins.filtered.Uniprot.fasta
$MAKER maker_functional_gff "$UNIPROT" blastp_uniprot.besthits.fullid.tsv "$GENES" \
    > filtered.genes.Uniprot.gff3

# 4. numbers for guiding questions -> blastp_summary.txt
{
    total=$(grep -c '>' "$PROTEINS")
    echo "Proteine total:        $total"
    for ref in uniprot TAIR10; do
        n=$(cut -f1 "blastp_${ref}.besthits.tsv" | sort -u | wc -l)
        awk -v n="$n" -v t="$total" -v r="$ref" 'BEGIN{printf "Mit Hit %-8s %6d (%.1f%%)\n", r":", n, 100*n/t}'
    done

    # how many of the best hits are "uncharacterized"?
    cut -f2 blastp_uniprot.besthits.fullid.tsv | sort -u > uniprot_hit_ids.txt
    unchar=$(grep '>' "$UNIPROT" | sed 's/^>//' | awk 'NR==FNR{h[$1];next} ($1 in h)' uniprot_hit_ids.txt - \
             | grep -ci 'uncharacterized' || true)
    echo "UniProt-Besthits mit 'Uncharacterized': $unchar"

    # length bias: proteins with vs. without UniProt hit
    echo; echo "Proteinlänge (aa): Gruppe  n  Median  Mittel"
    awk '/^>/{if(id)print id"\t"len; id=substr($1,2); len=0; next}{len+=length($0)} END{print id"\t"len}' "$PROTEINS" \
      | awk -F'\t' 'NR==FNR{hit[$1];next}{g=($1 in hit)?"hit":"no_hit"; print g"\t"$2}' \
            <(cut -f1 blastp_uniprot.besthits.tsv) - \
      | sort -k1,1 -k2,2n \
      | awk -F'\t' '{n[$1]++; v[$1,n[$1]]=$2; s[$1]+=$2}
            END{for(g in n){k=n[g]; m=(k%2)?v[g,(k+1)/2]:(v[g,k/2]+v[g,k/2+1])/2
                printf "  %-7s %6d %7d %7.0f\n", g, k, m, s[g]/k}}'

    echo; echo "FLC (AT5G10140) in TAIR10-Besthits:"
    grep -i 'AT5G10140' blastp_TAIR10.besthits.tsv || echo "  kein Besthit auf AT5G10140"
} > blastp_summary.txt
cat blastp_summary.txt