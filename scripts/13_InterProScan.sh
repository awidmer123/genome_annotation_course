#!/bin/bash
#SBATCH --job-name=InterPro
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/13_%x_%j.out
#SBATCH --error=logs/13_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

# setting paths
WORKDIR=/data/users/awidmer/genome_annotation_course
OUTDIR="$WORKDIR/output/09_MAKER"
FINALDIR="$OUTDIR/final"
COURSEDIR=/data/courses/genome_annotation_resources/2026
INTERPROSCAN_SIF="$COURSEDIR/containers/InterProScan_5.70.sif"
protein="assembly.all.maker.proteins.fasta"

# executing interproscan
apptainer exec --cleanenv --pwd "$FINALDIR" \
--bind $COURSEDIR/data/interproscan-5.70-102.0/data:/opt/interproscan/data \
--bind "$COURSEDIR" --bind "$WORKDIR" \
--bind "$FINALDIR/tmp":/temp \
"$INTERPROSCAN_SIF" \
/opt/interproscan/interproscan.sh \
-appl pfam --disable-precalc -f TSV \
--goterms --iprlookup --seqtype p \
-i ${protein}.renamed.fasta -o output.iprscan