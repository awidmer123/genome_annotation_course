#!/bin/bash
#SBATCH --job-name=Update
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/14_%x_%j.out
#SBATCH --error=logs/14_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

# setting paths
WORKDIR=/data/users/awidmer/genome_annotation_course
OUTDIR="$WORKDIR/output/09_MAKER"
FINALDIR="$OUTDIR/final"
COURSEDIR=/data/courses/genome_annotation_resources/2026
MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
gff="assembly.all.maker.noseq.gff"

cd "$FINALDIR"

# executing command
apptainer exec --cleanenv --env TMPDIR=/tmp --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$MAKER_SIF" ipr_update_gff ${gff}.renamed.gff output.iprscan > ${gff}.renamed.iprscan.gff