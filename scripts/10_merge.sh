#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=09:00:00
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --reservation=class-sbl-30004-2026-10-28
#SBATCH --output=logs/10_%x_%j.out
#SBATCH --error=logs/10_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

# set paths
COURSEDIR=/data/courses/genome_annotation_resources/2026
MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
OUTDIR="$WORKDIR/output/MAKER"

# create output directory
mkdir -p "$OUTDIR/tmp"
cd "$OUTDIR"

shopt -s nullglob
INDEXES=(*.maker.output/*_master_datastore_index.log)
test "${#INDEXES[@]}" -eq 1
INDEX="${INDEXES[0]}"


apptainer exec --cleanenv --env TMPDIR=/tmp --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" --
pwd "$OUTDIR" "$MAKER_SIF" \
gff3_merge -d "$INDEX" -o assembly.all.maker.gff


apptainer exec --cleanenv --env TMPDIR=/tmp --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" --
pwd "$OUTDIR" "$MAKER_SIF" \
gff3_merge -n -d "$INDEX" -o assembly.all.maker.noseq.gff


apptainer exec --cleanenv --env TMPDIR=/tmp --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" --
pwd "$OUTDIR" "$MAKER_SIF" \
fasta_merge -d "$INDEX" -o assembly