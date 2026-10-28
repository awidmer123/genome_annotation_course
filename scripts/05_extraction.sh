#!/bin/bash
#SBATCH --job-name=extraction
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --time=01:00:00
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --reservation=class-sbl-30004-2026-10-21
#SBATCH --output=logs/05_%x_%j.out
#SBATCH --error=logs/05_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch


# set paths
WORKDIR="/data/users/awidmer/genome_annotation_course"
COURSEDIR="/data/courses/genome_annotation_resources/2026"

# extracing copia and gipsy sequences
cd "$WORKDIR/output/01_EDTA"
genome="assembly_flye.fasta"

# creating output directory
mkdir -p "../05_extraction"

# Extract Copia sequences
apptainer exec -C --bind "$WORKDIR" "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
  seqkit grep -r -p "Copia" "$WORKDIR/output/01_EDTA/${genome}.mod.EDTA.TElib.fa" \
  > ../05_extraction/Copia_sequences.fa

# Extract Gypsy sequences
apptainer exec -C --bind "$WORKDIR" "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
  seqkit grep -r -p "Gypsy" "$WORKDIR/output/01_EDTA/${genome}.mod.EDTA.TElib.fa" \
  > ../05_extraction/Gypsy_sequences.fa