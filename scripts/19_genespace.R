library(GENESPACE)
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]   # GENESPACE-Verzeichnis mit bed/ und peptide/

gpar <- init_genespace(wd = wd, path2mcscanx = "/usr/local/bin", nCores = 20)
out <- run_genespace(gpar)

pangenome <- query_pangenes(out, bed = NULL, refGenome = "TAIR10", transform = TRUE,
                            showArrayMem = TRUE, showNSOrtho = TRUE, maxMem2Show = Inf)
saveRDS(pangenome, file = file.path(wd, "pangenome_matrix.rds"))
