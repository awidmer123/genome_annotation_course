#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --reservation=class-sbl-30004-2026-10-28
#SBATCH --output=logs/11_%x_%j.out
#SBATCH --error=logs/11_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

# setting paths
WORKDIR="/data/users/awidmer/genome_annotation_course"
OUTDIR="$WORKDIR/output/09_MAKER"
FINALDIR="$OUTDIR/final"
mkdir -p "$FINALDIR/tmp"
cd "$OUTDIR"

# specifying names
protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"
gff="assembly.all.maker.noseq.gff"

# copying
cp "$OUTDIR/$gff" "$FINALDIR/${gff}.renamed.gff"
cp "$OUTDIR/$protein" "$FINALDIR/${protein}.renamed.fasta"
cp "$OUTDIR/$transcript" "$FINALDIR/${transcript}.renamed.fasta"
cd "$FINALDIR"