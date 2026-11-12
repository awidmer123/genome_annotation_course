# 23_report_figures.R
# Builds the report figures that are not produced by the earlier pipeline steps
# and saves them (PDF + PNG) to report/graphics/.
#   fig02  TE superfamilies: % of genome (EDTA TEanno.sum)        + table tab02
#   fig03  Circos: Gypsy, Copia, TIR, Helitron and gene density
#   fig06  Cumulative AED distribution (+ Pfam share by AED class)
#   fig07  BUSCO: proteins vs. transcripts (longest isoform)
#   fig08  Protein length with vs. without UniProt hit
# Run from the project root (see 23_report_figures.sh).

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(readr); library(stringr); library(circlize)
})

# ---- paths -------------------------------------------------------------------
genome    <- "assembly_flye.fasta"
edta_dir  <- "output/01_EDTA"
final_dir <- "output/09_MAKER/final"
cmp_dir   <- "output/comparisons"
out_dir   <- "report/graphics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

te_sum   <- file.path(edta_dir, paste0(genome, ".mod.EDTA.TEanno.sum"))
te_gff   <- file.path(edta_dir, paste0(genome, ".mod.EDTA.TEanno.gff3"))
fai      <- file.path("data", paste0(genome, ".fai"))
gene_gff <- file.path(final_dir, "filtered.genes.renamed.gff3")
aed_txt  <- file.path(final_dir, "assembly.all.maker.renamed.gff.AED.txt")
pfam_txt <- file.path(final_dir, "pfam_by_AED.txt")
busco_p  <- Sys.glob(file.path(final_dir, "busco_proteins", "short_summary*.txt"))[1]
busco_t  <- Sys.glob(file.path(final_dir, "busco_transcripts", "short_summary*.txt"))[1]
proteins <- file.path(final_dir, "assembly.all.maker.proteins.fasta.renamed.filtered.fasta")
besthits <- file.path(cmp_dir, "blastp_uniprot.besthits.tsv")

theme_set(theme_bw(base_size = 12) +
          theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
save_both <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  message("saved ", name)
}
# colours per TE order, reused in fig02 and fig03
order_cols <- c(LTR = "#1b7837", TIR = "#2166ac", nonTIR = "#b2182b",
                LINE = "#e08214", Other = "grey50")

# ---- fig02: TE superfamilies -------------------------------------------------
sum_lines <- readLines(te_sum)
sum_lines <- sum_lines[seq_len(grep("^Repeat Stats", sum_lines)[1] - 1)]  # first block only
genome_bp <- as.numeric(str_match(sum_lines[grep("Total Length", sum_lines)[1]], "([0-9]+) bp")[, 2])
te <- tibble(line = sum_lines) |>
  mutate(order = ifelse(str_detect(line, "^[A-Za-z_]+\\s+--"), str_extract(line, "^[A-Za-z_]+"), NA)) |>
  fill(order) |>
  filter(str_detect(line, "^\\s*\\S+\\s+[0-9]+\\s+[0-9]+\\s+[0-9.]+%"),
         !str_detect(line, "^\\s*(Total|total interspersed)")) |>
  mutate(fields = str_split(str_trim(line), "\\s+"),
         superfamily = sapply(fields, `[`, 1),
         count       = as.integer(sapply(fields, `[`, 2)),
         bp          = as.numeric(sapply(fields, `[`, 3)),
         pct         = 100 * bp / genome_bp) |>
  select(order, superfamily, count, bp, pct)
# repeat_fragment / snRNA are top-level lines without an order header
te <- te |> mutate(order = ifelse(superfamily %in% c("repeat_fragment", "snRNA"), "Other", order))
write_csv(te |> mutate(pct = round(pct, 2)), file.path(out_dir, "tab02_TE_superfamilies.csv"))

te_plot <- te |> filter(superfamily != "snRNA") |>
  mutate(label = ifelse(superfamily == "unknown", "LTR unknown", superfamily),
         label = reorder(label, pct),
         order = factor(order, levels = names(order_cols)))
p2 <- ggplot(te_plot, aes(pct, label, fill = order)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%.2f%%", pct)), hjust = -0.15, size = 3.4) +
  scale_fill_manual(values = order_cols, name = "Order") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "% of assembly (bp masked)", y = NULL,
       title = sprintf("TE content of Mh-0: %.1f%% of %.1f Mb",
                       sum(te$bp[te$superfamily != "snRNA"]) / genome_bp * 100, genome_bp / 1e6))
save_both(p2, "fig02_TE_superfamilies", 8, 5)

# ---- fig03: circos with TEs and genes ----------------------------------------
fa_idx <- read_tsv(fai, col_names = FALSE, show_col_types = FALSE) |> select(chr = X1, len = X2)
keep <- fa_idx |> arrange(desc(len)) |> filter(len >= 1e6)        # scaffolds >= 1 Mb
if (nrow(keep) < 5) keep <- fa_idx |> arrange(desc(len)) |> slice_head(n = 10)
ideo <- data.frame(chr = keep$chr, start = 0, end = keep$len)

read_gff <- function(f) read_tsv(f, comment = "#", col_names = FALSE, show_col_types = FALSE,
                                 col_types = cols(.default = "c")) |>
  transmute(chr = X1, type = X3, start = as.numeric(X4), end = as.numeric(X5)) |>
  filter(chr %in% keep$chr)
tes   <- read_gff(te_gff)
genes <- read_gff(gene_gff) |> filter(type == "gene")
bed_of <- function(d) as.data.frame(d[, c("chr", "start", "end")])
tracks <- list(
  Gypsy    = list(d = filter(tes, type == "Gypsy_LTR_retrotransposon"), col = "#1b7837"),
  Copia    = list(d = filter(tes, type == "Copia_LTR_retrotransposon"), col = "#7fbc41"),
  TIR      = list(d = filter(tes, str_detect(type, "TIR_transposon")),  col = "#2166ac"),
  Helitron = list(d = filter(tes, type == "helitron"),                  col = "#b2182b"),
  Genes    = list(d = genes,                                            col = "black"))

draw_circos <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.after = 1, track.margin = c(0.005, 0.005),
             cell.padding = c(0, 0, 0, 0))
  circos.genomicInitialize(ideo, plotType = c("axis", "labels"), axis.labels.cex = 0.45,
                           labels.cex = 0.6, major.by = 5e6)
  for (nm in names(tracks)) {
    circos.genomicDensity(bed_of(tracks[[nm]]$d), count_by = "number", window.size = 1e5,
                          col = tracks[[nm]]$col, track.height = 0.1)
  }
  legend("bottomleft", inset = 0.01, bty = "n", cex = 0.8, title = "outer to inner",
         legend = names(tracks), fill = sapply(tracks, `[[`, "col"))
  circos.clear()
}
pdf(file.path(out_dir, "fig03_circos_TE_genes.pdf"), width = 9, height = 9); draw_circos(); dev.off()
png(file.path(out_dir, "fig03_circos_TE_genes.png"), width = 9, height = 9, units = "in", res = 300)
draw_circos(); dev.off()
message("saved fig03_circos_TE_genes (", nrow(keep), " scaffolds >= 1 Mb shown)")

# ---- fig06: cumulative AED ----------------------------------------------------
aed <- read_tsv(aed_txt, show_col_types = FALSE)
names(aed) <- c("AED", "cum_fraction")
pfam <- readLines(pfam_txt)
p6 <- ggplot(aed, aes(AED, cum_fraction)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_step(linewidth = 1, colour = "#2166ac") +
  annotate("label", x = 0.98, y = 0.08, hjust = 1, vjust = 0, size = 3.3, label.size = 0,
           label = paste(c("Pfam domain per AED class:", gsub("mit Pfam", "with Pfam", pfam)), collapse = "\n")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Annotation Edit Distance (AED)", y = "Cumulative fraction of mRNAs",
       title = "Evidence support of MAKER gene models")
save_both(p6, "fig06_AED_cdf", 7, 5)

# ---- fig07: BUSCO --------------------------------------------------------------
read_busco <- function(f, label) {
  s <- str_match(paste(readLines(f), collapse = " "),
                 "C:([0-9.]+)%\\[S:([0-9.]+)%,D:([0-9.]+)%\\],F:([0-9.]+)%,M:([0-9.]+)%,n:([0-9]+)")
  tibble(set = label, n = as.integer(s[7]),
         category = c("Complete single-copy", "Complete duplicated", "Fragmented", "Missing"),
         pct = as.numeric(s[3:6]))
}
busco <- bind_rows(read_busco(busco_p, "Proteins"), read_busco(busco_t, "Transcripts")) |>
  mutate(category = factor(category, levels = c("Missing", "Fragmented",
                                                "Complete duplicated", "Complete single-copy")))
busco_lab <- busco |> group_by(set) |>
  summarise(lab = sprintf("C:%.1f%% [S:%.1f%%, D:%.1f%%], F:%.1f%%, M:%.1f%%, n=%d",
                          sum(pct[str_detect(category, "Complete")]),
                          pct[category == "Complete single-copy"], pct[category == "Complete duplicated"],
                          pct[category == "Fragmented"], pct[category == "Missing"], first(n)))
p7 <- ggplot(busco, aes(pct, set, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(data = busco_lab, aes(x = 2, y = set, label = lab), inherit.aes = FALSE,
            hjust = 0, colour = "white", size = 3.4, fontface = "bold") +
  scale_fill_manual(values = c("Complete single-copy" = "#2166ac", "Complete duplicated" = "#67a9cf",
                               "Fragmented" = "#fdb863", "Missing" = "#d6604d"), name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "% of BUSCO groups (brassicales_odb10)", y = NULL,
       title = "Completeness of the Mh-0 annotation (longest isoform per gene)") +
  theme(legend.position = "bottom")
save_both(p7, "fig07_BUSCO", 9, 3.2)

# ---- fig08: protein length with vs. without UniProt hit -----------------------
fa <- readLines(proteins)
hdr <- grepl("^>", fa)
ids <- sub("^>(\\S+).*", "\\1", fa[hdr])
lens <- tapply(nchar(fa[!hdr]), cumsum(hdr)[!hdr], sum)
prot <- tibble(id = ids, length = as.numeric(lens))
hit_ids <- unique(read_tsv(besthits, col_names = FALSE, show_col_types = FALSE)$X1)
prot <- prot |> mutate(group = ifelse(id %in% hit_ids, "UniProt hit", "No UniProt hit"))
stats <- prot |> group_by(group) |>
  summarise(n = n(), median = median(length), .groups = "drop") |>
  mutate(lab = sprintf("n = %s\nmedian = %d aa", format(n, big.mark = "'"), as.integer(median)))
p8 <- ggplot(prot, aes(group, length, fill = group)) +
  geom_violin(scale = "width", alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_text(data = stats, aes(y = 4500, label = lab), size = 3.4, vjust = 1) +
  scale_y_log10(breaks = c(30, 100, 300, 1000, 3000)) +
  scale_fill_manual(values = c("UniProt hit" = "#2166ac", "No UniProt hit" = "#d6604d"), guide = "none") +
  labs(x = NULL, y = "Protein length (aa, log scale)",
       title = sprintf("%.1f%% of %s proteins have a UniProt hit",
                       100 * mean(prot$group == "UniProt hit"), format(nrow(prot), big.mark = "'")))
save_both(p8, "fig08_protein_length_by_hit", 6, 5)

message("done – figures in ", out_dir)