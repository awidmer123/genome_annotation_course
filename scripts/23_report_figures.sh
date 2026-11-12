#!/bin/bash
#SBATCH --job-name=report_figs
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/23_%x_%j.out
#SBATCH --error=logs/23_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch

# Builds the report figures (fig02, 03, 06, 07, 08) into report/graphics/.
# Submit from the project root: sbatch scripts/23_report_figures.sh
set -euo pipefail

WORKDIR=/data/users/$USER/genome_annotation_course
COURSEDIR=/data/courses/genome_annotation_resources/2026
cd "$WORKDIR"

apptainer exec --cleanenv --bind "$WORKDIR" --bind "$COURSEDIR" --pwd "$WORKDIR" \
    "$COURSEDIR/containers/Annotation_analysis_R4.4.sif" \
    Rscript --vanilla scripts/23_report_figures.R