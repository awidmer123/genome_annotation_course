#!/bin/bash
#SBATCH --job-name=TE_sorter_plots
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --reservation=class-sbl-30004-2026-10-14
#SBATCH --output=logs/03_%x_%j.out
#SBATCH --error=logs/03_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

# set paths
WORKDIR="/data/users/awidmer/genome_annotation_course"
COURSEDIR="/data/courses/genome_annotation_resources/2026"

# execute the prepared R script
mkdir -p plots
cd "$WORKDIR/output/02_LTR_full_length"
apptainer exec -C --bind "$WORKDIR" --bind "$COURSEDIR" --pwd "$WORKDIR/output/02_LTR_full_length" \
"$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
Rscript --vanilla "$WORKDIR/scripts/03_LTR_identity.R"
