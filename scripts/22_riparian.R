library(GENESPACE)
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
focal_genome <- args[2]
# Scaffolds to mirror: ab
flip_chrs <- if (length(args) > 2) args[-(1:2)] else character(0)

load(file.path(wd, "results", "gsParams.rda"))   # lädt gsParam
outdir <- file.path(wd, "processed_riparian")
dir.create(outdir, showWarnings = FALSE)

flips <- if (length(flip_chrs)) data.frame(genome = focal_genome, chr = flip_chrs) else NULL

# 1. Mh-0 vs. TAIR10
plot_riparian(gsParam, genomeIDs = c("TAIR10", focal_genome), refGenome = "TAIR10",
              reorderBySynteny = TRUE,
              pdfFile = file.path(outdir, "my_accession_vs_TAIR10.pdf"))

# 2. With mirrored scaffolds
plot_riparian(gsParam, genomeIDs = c("TAIR10", focal_genome), refGenome = "TAIR10",
              reorderBySynteny = TRUE, syntenyWeight = 1, invertTheseChrs = flips,
              pdfFile = file.path(outdir, "my_accession_vs_TAIR10_oriented.pdf"))

# 3. All five genomes
plot_riparian(gsParam, refGenome = "TAIR10", reorderBySynteny = TRUE, syntenyWeight = 1,
              invertTheseChrs = flips, scalePlotWidth = 1.5,
              pdfFile = file.path(outdir, "all_genomes_TAIR10.pdf"))
