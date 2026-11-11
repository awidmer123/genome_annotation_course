#!/bin/bash
#SBATCH --job-name=consistent_ids
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/12_%x_%j.out
#SBATCH --error=logs/12_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

# set paths
WORKDIR="/data/users/awidmer/genome_annotation_course"
OUTDIR="$WORKDIR/output/09_MAKER"
COURSEDIR=/data/courses/genome_annotation_resources/2026
FINALDIR="$OUTDIR/final"
MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
prefix="ALT"

# setting names
protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"
gff="assembly.all.maker.noseq.gff"

cd "$FINALDIR"

apptainer exec --cleanenv --env TMPDIR=/tmp --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$MAKER_SIF" maker_map_ids --prefix $prefix --justify 7 ${gff}.renamed.gff > id.map

apptainer exec --cleanenv --env TMPDIR=/tmp --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$MAKER_SIF" map_gff_ids id.map ${gff}.renamed.gff

apptainer exec --cleanenv --env TMPDIR=/tmp --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$MAKER_SIF" map_fasta_ids id.map ${protein}.renamed.fasta

apptainer exec --cleanenv --env TMPDIR=/tmp --bind "$OUTDIR/tmp":/tmp --bind "$COURSEDIR" --bind "$WORKDIR" --pwd "$FINALDIR" "$MAKER_SIF" map_fasta_ids id.map ${transcript}.renamed.fasta