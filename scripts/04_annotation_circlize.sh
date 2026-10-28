#!/bin/bash
#SBATCH --job-name=circlize
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --reservation=class-sbl-30004-2026-10-21
#SBATCH --output=logs/04_%x_%j.out
#SBATCH --error=logs/04_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

# IMPORTANT NOTE:
# first run the following command in an interactive session to guarantee a successful execution:
# apptainer exec -C --bind "$WORKDIR" "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" samtools faidx "$WORKDIR"/data/assembly_flye.fasta

mkdir -p "output/04_circlize"

# set paths
WORKDIR="/data/users/awidmer/genome_annotation_course"
COURSEDIR="/data/courses/genome_annotation_resources/2026"

# execute the prepared R script
apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR"
"$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
Rscript --vanilla "$WORKDIR/scripts/04_annotation_circlize.R"
