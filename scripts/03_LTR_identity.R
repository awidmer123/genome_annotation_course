# In this script, we analyze the identity of full-length LTR retrotransposons (LTR-RTs)
# annotated in a genome. We read a GFF file with LTR-RT annotations,
# extract the LTR identity values, merge with classification data from TEsorter,
# and generate plots showing the distribution of LTR identities per clade within
# the Copia and Gypsy superfamilies.

library(tidyverse)
library(data.table)
library(cowplot)

#-------------------------------------------------
# Input files (edit to match your data raw LTR GFF and TEsorter classification output)
#-------------------------------------------------
gff_file <- "intact_LTR.gff3"
cls_file <- "intact_LTR.fa.rexdb-plant.cls.tsv"
# Run this script in output/LTR_full_length after TEsorter on intact LTR elements.
#-------------------------------------------------
# Read and preprocess input data
#-------------------------------------------------
message("Reading GFF: ", gff_file)
anno <- read.table(gff_file, sep = "\t", header = FALSE)

# Remove subfeatures (Terminal repeats, TSDs) so we keep top-level TE annotations
exclude_feats <- c("long_terminal_repeat", "repeat_region", "target_site_duplication")
anno <- anno %>% filter(!V3 %in% exclude_feats)

head(anno)

# Extract Name and ltr_identity from the ninth column (attributes). This uses regex.
anno <- anno %>%
  mutate(
    # extract Name=... from attributes (V9)
    Name = str_match(V9, "Name=([^;]+)")[, 2],
    # extract ltr_identity=...
    Identity = as.numeric(str_match(V9, "ltr_identity=([^;]+)")[, 2]),
    # compute length using inclusive GFF coordinates
    length = as.numeric(V5) - as.numeric(V4) + 1
  ) %>%
  # keep only the columns used downstream
  select(V1, V4, V5, V3, Name, Identity, length)

head(anno)

message("Reading classification: ", cls_file)
# Read classification table (TE name in first column). If your file doesn't have a header, set header=FALSE.
cls <- fread(cls_file, sep = "\t", header = TRUE)
setnames(cls, 1, "TE")

head(cls)

# Intact LTR headers use element-coordinate IDs followed by #Classification.
# Split on '#'; the first part is the Name, the second part is the Classification.
cls <- cls %>%
  separate(TE, into = c("Name", "Classification"), sep = "#", fill = "right")

head(cls)

## Merge annotation with classification table
# Use a left join so all annotated TEs are kept even if they have no classification match
anno_cls <- merge(anno, cls, by = "Name", all.x = TRUE)

head(anno_cls)

# Quick checks: how many per Superfamily/Clade (may be NA if classification missing)
message("Counts per Superfamily")
print(table(anno_cls$Superfamily, useNA = "ifany"))
message("Counts per Clade")
print(table(anno_cls$Clade, useNA = "ifany"))


#-------------------------------------------------
# Plot setup
#-------------------------------------------------
# binwidth controls histogram resolution around identity values
binwidth <- 0.005
# x axis limits, these can be adjusted if needed, minimum identity in your data may differ
xlims <- c(0.80, 1.00)
bin_breaks <- seq(xlims[1], xlims[2], by = binwidth)
sf_colors <- c(Copia = "#1b9e77", Gypsy = "#d95f02")
plot_data <- anno_cls %>%
  filter(Superfamily %in% names(sf_colors), !is.na(Identity))

# Compute a single y-max across ALL Copia and Gypsy clades. This ensures consistent y-axis scaling.
global_ymax <- plot_data %>%
  filter(!is.na(Clade), between(Identity, xlims[1], xlims[2])) %>%
  # bin Identity into consistent breaks and count occurrences
  count(Superfamily, Clade, Identity = cut(Identity, bin_breaks, include.lowest = TRUE)) %>%
  pull(n) %>%
  max(na.rm = TRUE)


#-------------------------------------------------
# Plot function for one superfamily
#-------------------------------------------------
plot_by_clade <- function(df, sf, ymax) {
  df %>%
    filter(Superfamily == sf, !is.na(Identity), !is.na(Clade)) %>%
    ggplot(aes(x = Identity)) +
    # histogram with color coding per superfamily
    geom_histogram(
      breaks = bin_breaks,
      color = "black",
      fill = sf_colors[[sf]]
    ) +
    # One panel per clade within this superfamily, stacked in one column.
    # All panels share the same axis scales by default.
    facet_wrap(~Clade, ncol = 1) +
    # x axis focused around xlims values
    scale_x_continuous(breaks = seq(xlims[1], xlims[2], 0.05)) +
    # set y limit to the provided ymax
    coord_cartesian(xlim = xlims, ylim = c(0, ymax), expand = FALSE) +
    theme_cowplot() +
    theme(
      strip.background = element_rect(fill = "#f0f0f0"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = sf, x = "Identity", y = "Count")
}

#-------------------------------------------------
# Generate Copia and Gypsy plots
#-------------------------------------------------
p_copia <- plot_by_clade(plot_data, "Copia", global_ymax)
p_gypsy <- plot_by_clade(plot_data, "Gypsy", global_ymax)

# Combine with cowplot side-by-side
combined <- plot_grid(p_copia, p_gypsy, ncol = 2)

out_dir <- "../../plots/02_LTR"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave("../../plots/02_LTR/01_LTR_Copia_Gypsy_cladelevel.png", combined, width = 12, height = 10, dpi = 300)
ggsave("../../plots/02_LTR/01_LTR_Copia_Gypsy_cladelevel.pdf", combined, width = 12, height = 10)

#-------------------------------------------------
# Global Copia versus Gypsy comparison (all clades pooled)
#-------------------------------------------------
global_plot <- plot_data %>%
  ggplot(aes(x = Identity, fill = Superfamily)) +
  geom_histogram(breaks = bin_breaks, color = "black") +
  # Pool all clades into two panels: Copia and Gypsy, side by side.
  # Both panels share the same axis scales by default.
  facet_wrap(~Superfamily, ncol = 2) +
  scale_fill_manual(values = sf_colors) +
  scale_x_continuous(breaks = seq(xlims[1], xlims[2], 0.05)) +
  coord_cartesian(xlim = xlims, expand = FALSE) +
  theme_cowplot() +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    legend.position = "none"
  ) +
  labs(x = "Identity", y = "Count", title = "Full-length LTRs: all clades pooled")

ggsave("../../plots/02_LTR/02_LTR_Copia_Gypsy_global.png", global_plot, width = 10, height = 5, dpi = 300)
ggsave("../../plots/02_LTR/02_LTR_Copia_Gypsy_global.pdf", global_plot, width = 10, height = 5)
