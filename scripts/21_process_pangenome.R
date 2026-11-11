library(tidyverse)


# ============================================================================
# Settings for your GENESPACE run
# ============================================================================
wd <- "/data/users/awidmer/genome_annotation_course/output/comparisons/GENESPACE"
focal_genome <- "Mh_0"

options(tibble.width = Inf)

# Keep all results from this script together inside the GENESPACE directory.
output_dir <- file.path(wd, "process_pangenome")
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# STEP 1: Read the TAIR10-reference pan-gene table.
# ============================================================================
# pgID = positional pan-gene entry; genome/id = accession and gene identifier.
# flag: PASS marks syntenic genes; array marks additional tandem-array members.
# The syntenic set includes both PASS and array; non-syntenic matches are excluded.
# chr/start/end are gene coordinates; interpChr/interpOrd are TAIR10-relative.
pg <- read_tsv(
  file.path(wd, "pangenes/TAIR10_pangenes.txt.gz"),
  show_col_types = FALSE
) %>% filter(!is.na(interpOrd)) # removes records without reference-relative positions.
message("\nFirst six rows of the pangenome table (entries with a reference-relative position).")
print(head(pg))

genomes <- sort(unique(pg$genome))
n_genomes <- length(genomes)
message("\nGenomes included in the pangenome table:")
print(genomes)

genes <- pg %>%
  filter(!is.na(id), nzchar(id)) %>%
  distinct(pgID, genome, id)

# ============================================================================
# STEP 2: How many copies does each entry have in each accession?
# ============================================================================
# Count the recorded genes in each entry/accession.
counts <- genes %>%
  count(pgID, genome, name = "copies")

# List every entry/accession combination, including those with no matches.
all_combinations <- expand_grid(
  pgID = unique(pg$pgID),
  genome = genomes
)

# Attach the counts and replace missing counts with zero.
counts <- all_combinations %>%
  left_join(counts, by = c("pgID", "genome")) %>%
  mutate(copies = replace_na(copies, 0L))
message("\nFirst six entry/accession copy counts.")
print(head(counts))

# ============================================================================
# STEP 3: Summarise presence and plot its frequency.
# ============================================================================
# Core = all accessions; accession-specific = one; accessory = some but not all.
entries <- counts %>%
  group_by(pgID) %>%
  summarise(
    n_present = sum(copies > 0),
    single_copy_core = all(copies == 1),
    present_in_focal = any(genome == focal_genome & copies > 0),
    .groups = "drop"
  ) %>%
  mutate(category = case_when(
    n_present == 0 ~ "no_recorded_match",
    n_present == n_genomes ~ "core",
    n_present == 1 ~ "accession_specific",
    TRUE ~ "accessory"
  ))
message("\nNumber of pan-gene entries in each presence category.")
print(entries %>% count(category))

# A gene can match several entries. Count distinct IDs for unique gene totals.

genome_summary <- genes %>%
  distinct(genome, id) %>%
  count(genome, name = "unique_genes")

# Include every genome, even if it has no valid gene IDs in this table.
genome_summary <- expand_grid(genome = genomes) %>%
  left_join(genome_summary, by = "genome") %>%
  mutate(unique_genes = replace_na(unique_genes, 0L))
message("\nUnique gene IDs represented in the selected pangenome matches, per accession.")
print(genome_summary)

frequency <- entries %>% count(n_present, name = "pan_gene_entries")
frequency_plot <- ggplot(frequency, aes(n_present, pan_gene_entries)) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(breaks = 0:n_genomes) +
  labs(
    x = "Number of accessions", y = "Pan-gene entries",
    title = "Pangenome presence frequency",
    subtitle = "All recorded matches"
  ) +
  theme_minimal()
message("\nPlot: number of pan-gene entries represented in each number of accessions.")
print(frequency_plot)

# Save the PDF in process_pangenome/plots.
ggsave(file.path(plot_dir, "pangenome_presence_frequency.pdf"),
  plot = frequency_plot, width = 7, height = 5
)


# ============================================================================
# STEP 4: Count PASS and array copies, including zeros for missing matches.
# ============================================================================
syntenic_counts <- pg %>%
  filter(flag %in% c("PASS", "array"), !is.na(id), nzchar(id)) %>%
  distinct(pgID, genome, id) %>%
  count(pgID, genome, name = "copies")

syntenic_counts <- all_combinations %>%
  left_join(syntenic_counts, by = c("pgID", "genome")) %>%
  mutate(copies = replace_na(copies, 0L))

# ============================================================================
# STEP 5: Identify copy-variable core entries separately for each gene set.
# ============================================================================
comparison_counts <- bind_rows(
  counts %>% mutate(gene_set = "All recorded genes"),
  syntenic_counts %>% mutate(gene_set = "Syntenic + array genes (PASS + array)")
)

copy_variable_core <- comparison_counts %>%
  group_by(gene_set, pgID) %>%
  summarise(
    core = all(copies > 0),
    min_copies = min(copies),
    max_copies = max(copies),
    .groups = "drop"
  ) %>%
  filter(core, max_copies > min_copies)

# Report the number of qualifying entries in each set, including zero.
core_summary <- tibble(gene_set = unique(comparison_counts$gene_set)) %>%
  left_join(count(copy_variable_core, gene_set, name = "entries"),
    by = "gene_set"
  ) %>%
  mutate(entries = replace_na(entries, 0L))
message("\nNumber of copy-variable core entries in each gene set (all recorded genes versus PASS + array genes).")
print(core_summary)

# ============================================================================
# STEP 6: Plot copy-number differences for up to 20 entries per set.
# ============================================================================
# Rank each gene set separately by copies in the focal accession.
# Ties are resolved by pgID; the two panels may contain different entries.
focal_counts <- comparison_counts %>%
  filter(genome == focal_genome) %>%
  select(gene_set, pgID, focal_copies = copies)

plot_entries <- copy_variable_core %>%
  left_join(focal_counts, by = c("gene_set", "pgID")) %>%
  group_by(gene_set) %>%
  arrange(desc(focal_copies), pgID, .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(focal_rank = row_number()) %>%
  ungroup()

# Separate row labels allow a different pgID order in each panel.
plot_order <- paste(plot_entries$gene_set, plot_entries$pgID, sep = "::")
plot_counts <- comparison_counts %>%
  semi_join(plot_entries, by = c("gene_set", "pgID")) %>%
  mutate(
    plot_row = factor(paste(gene_set, pgID, sep = "::"),
      levels = rev(plot_order)
    ),
    gene_set = factor(gene_set, levels = core_summary$gene_set)
  )

if (nrow(plot_counts) > 0) {
  cnv_plot <- ggplot(plot_counts, aes(genome, plot_row, fill = copies)) +
    geom_tile(colour = "white") +
    # White numbers on dark tiles; black numbers on the brighter tiles.
    geom_text(
      aes(
        label = copies,
        colour = copies > min(copies) + 0.6 * (max(copies) - min(copies))
      ),
      size = 4, fontface = "bold"
    ) +
    scale_colour_manual(
      values = c("FALSE" = "white", "TRUE" = "black"),
      guide = "none"
    ) +
    facet_wrap(~gene_set, scales = "free_y", drop = FALSE) +
    scale_y_discrete(labels = function(x) sub("^.*::", "", x)) +
    scale_fill_viridis_c() +
    labs(
      x = "Accession", y = "Pan-gene entry (pgID)", fill = "Copies",
      title = "Copy-variable core entries",
      subtitle = paste("Top 20 per set, ranked by copies in", focal_genome, "(highest at top)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  message("\nPlot: gene copy numbers across accessions for the top 20 copy-variable core entries ranked by focal-accession copies in each gene set.")
  print(cnv_plot)
  ggsave(file.path(plot_dir, "copy_variable_core_comparison.pdf"),
    cnv_plot,
    width = 11, height = 9
  )
} else {
  message("No copy-variable core entries were found in either set.")
}

# Extra copies outside the conserved position could be genuine duplications,
# but also assembly artifacts.
# To investigate a few high-copy entries, filter the pg table by pgID,
# select the TAIR10 rows, and check their chr values:
# ChrC or ChrM indicate chloroplast or mitochondrial reference genes.


# Save tables for inspecting individual entries.

# List the actual genes in your accession and TAIR10, including all flags.
# Each row describes one gene's association with a pgID.
gene_positions <- pg %>%
  filter(genome %in% c(focal_genome, "TAIR10"), !is.na(id), nzchar(id)) %>%
  select(pgID, genome, id, flag, chr, start, end) %>%
  distinct() %>%
  arrange(pgID, genome, chr, start)
write_tsv(gene_positions, file.path(table_dir, "focal_TAIR10_genes.tsv"))
