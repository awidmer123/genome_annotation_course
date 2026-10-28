#!/bin/bash
#SBATCH --job-name=Maker_gene_annotation
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=09:00:00
#SBATCH --account=class-sbl-30004-ws-2026
#SBATCH --reservation=class-sbl-30004-2026-10-28
#SBATCH --output=logs/08_%x_%j.out
#SBATCH --error=logs/08_%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=andri.widmer@unifr.ch


set -euo pipefail

COURSEDIR=/data/courses/genome_annotation_resources/2026
WORKDIR="$SLURM_SUBMIT_DIR"
OUTDIR="$WORKDIR/output/MAKER"
MAKER_SIF="$COURSEDIR/containers/MAKER_3.01.03.sif"
AUGUSTUS_CONFIG_PATH="$COURSEDIR/tools/MAKER/augustus_config"

# Before submitting, copy the four templates from:
# /data/courses/genome_annotation_resources/2026/tools/MAKER/control_files/
# into your project root, then edit maker_opts.ctl:
# genome = your genome FASTA; est = your Trinity FASTA; rmlib = your EDTA TE library.
# Protein references already point to the course directory.
# Submit from the project root: sbatch scripts/XX-run_maker.sh
# The copy below uses YOUR EDITED project files, then runs MAKER in output/MAKER.
unset SINGULARITY_BIND APPTAINER_BIND SINGULARITY_TMPDIR
mkdir -p "$OUTDIR/tmp"
cp "$WORKDIR/maker_opts.ctl" "$WORKDIR/maker_bopts.ctl" \
    "$WORKDIR/maker_exe.ctl" "$WORKDIR/maker_evm.ctl" "$OUTDIR/"
cd "$OUTDIR"



# Use the container's MPICH launcher; keep cpus=1 in maker_opts.ctl.
apptainer exec --cleanenv --env AUGUSTUS_CONFIG_PATH="$AUGUSTUS_CONFIG_PATH" \
    --bind /data --bind "$WORKDIR" --bind "$OUTDIR/tmp":/TMP \
    --pwd "$OUTDIR" "$MAKER_SIF" \
    mpiexec -launcher fork -n "$SLURM_CPUS_PER_TASK" \
    maker -mpi --ignore_nfs_tmp -TMP /TMP \
    maker_opts.ctl maker_bopts.ctl maker_exe.ctl maker_evm.ctl



