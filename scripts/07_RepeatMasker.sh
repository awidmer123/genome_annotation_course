#!/bin/bash
#SBATCH --job-name=parseRM
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --time=01:00:00
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --reservation=class-sbl-30004-2026-10-21
#SBATCH --output=logs/07_%x_%j.out
#SBATCH --error=logs/07_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

set -euo pipefail

# genome name
genome="assembly_flye.fasta"

# base paths (zuerst!)
WORKDIR="/data/users/awidmer/genome_annotation_course"
COURSEDIR="/data/courses/genome_annotation_resources/2026"

# derived paths
ANNODIR="$WORKDIR/output/01_EDTA/${genome}.mod.EDTA.anno"
OUTDIR="$WORKDIR/output/07_parseRM"
RM_OUT="$ANNODIR/${genome}.mod.EDTA.RM.out"
RM_NONZERO="$OUTDIR/${genome}.mod.EDTA.RM.nonzero.out"
CONTAINER="$COURSEDIR/containers/Annotation_analysis_R4.4.sif"

mkdir -p "$OUTDIR"

# remove 0% divergence hits (keep header)
awk 'NR == 1 || $2 != 0' "$RM_OUT" > "$RM_NONZERO"

# bin divergence (output lands next to $RM_NONZERO, i.e. in $OUTDIR)
apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR" "$CONTAINER" \
    perl "$WORKDIR/scripts/05_parseRM.pl" -i "$RM_NONZERO" -l 50,1 -v

