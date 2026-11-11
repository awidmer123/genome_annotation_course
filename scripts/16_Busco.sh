#!/bin/bash
#SBATCH --job-name=Busco
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/16_%x_%j.out
#SBATCH --error=logs/16_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

set -euo pipefail

# --- same variables as in the filter script ---
WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
OUTDIR=$WORKDIR/output/09_MAKER
FINALDIR=$OUTDIR/final
ANALYSIS_SIF=$COURSEDIR/containers/Annotation_analysis_R4.4.sif
AGAT_SIF=$COURSEDIR/containers/AGAT_v1.4.1.sif

protein=assembly.all.maker.proteins.fasta
transcript=assembly.all.maker.transcripts.fasta
# ----------------------------------------------

cd "$FINALDIR"
ANALYSIS="apptainer exec --cleanenv --bind $COURSEDIR --bind $WORKDIR --pwd $FINALDIR $ANALYSIS_SIF"

# 0. Longest isoform per gene (gene = everything before -R*, ties: first one wins)
pick_longest () {
    local in=$1 out=$2
    $ANALYSIS seqkit fx2tab -n -i -l "$in" \
      | awk -F'\t' '{ id=$1; g=id; sub(/-R[A-Z]+$/, "", g)
                      if (!(g in len) || $2 > len[g]) { len[g]=$2; best[g]=id } }
                    END { for (g in best) print best[g] }' > "${out}.ids"
    $ANALYSIS seqkit grep -f "${out}.ids" "$in" -o "$out"
    echo "$out: $(grep -c '>' "$out") sequences (from $(grep -c '>' "$in"))"
}
pick_longest ${protein}.renamed.filtered.fasta    ${protein}.renamed.longest.fasta
pick_longest ${transcript}.renamed.filtered.fasta ${transcript}.renamed.longest.fasta

# 1. BUSCO
module load BUSCO/5.4.2-foss-2021a
export PYTHONNOUSERSITE=1
busco -i ${protein}.renamed.longest.fasta    -l brassicales_odb10 -o busco_proteins    -m proteins      -c 4 -f
busco -i ${transcript}.renamed.longest.fasta -l brassicales_odb10 -o busco_transcripts -m transcriptome -c 4 -f

# 3. AGAT statistics
apptainer exec --cleanenv --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$AGAT_SIF" \
    agat_sp_statistics.pl --gff filtered.genes.renamed.gff3 --output annotation.stat