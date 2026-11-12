# Genome annotation of *Arabidopsis thaliana* accession Mh-0

Follow-up to the genome assembly course ([GenomeTranscriptomeAssembly](https://github.com/awidmer123/GenomeTranscriptomeAssembly)).
Starting from the Flye assembly of PacBio HiFi reads of accession **Mh-0** (Mühlen, Poland), this repository annotates
transposable elements (TEs) and protein-coding genes, assigns putative functions, and compares Mh-0 with TAIR10 (Col-0)
and three other accessions (Taz-0, Altai-5, Ice-1; Lian et al.) in a pangenome / synteny framework.

All jobs were run on the UNIBE cluster of the bioinformatic unit using the course containers in
`/data/courses/genome_annotation_resources/2026`.

## Directory layout

```
data/      assembly_flye.fasta (input, not tracked)
scripts/   numbered scripts, run in order (see below)
output/    tool outputs, one folder per step (not tracked)
plots/     figures (not tracked)
logs/      SLURM logs (not tracked)
maker_*.ctl  MAKER control files (paths to genome, Trinity transcripts, EDTA TE library, protein evidence)
```

## Running order

Submit every `.sh` script from the project root: `sbatch scripts/<script>.sh`.
`.R` scripts are called by the `.sh` script with the same number.

| Step | Scripts | What it does | Main output |
|---|---|---|---|
| **1. TE annotation** | `01_edta.sh` | EDTA de novo TE library + genome-wide TE annotation (`--sensitive 1 --anno 1`, TAIR10 CDS to exclude genes) | `output/01_EDTA/*.TEanno.gff3`, `*.TEanno.sum`, `*.TElib.fa` |
| | `02_TEsorter.sh`, `03_LTR_identity.R` | Clade classification of intact LTR-RTs (TEsorter, REXdb-plant) and LTR identity (insertion age) per clade | `plots/02_LTR/` |
| | `04_annotation_circlize.sh/.R` | Circos plot of TE superfamily density along the scaffolds | `output/04_circlize/` |
| | `05_extraction.sh`, `06_TESorter._II.sh` | Extract Copia/Gypsy from the TE library and classify them into clades | `output/06_TESorter/` |
| | `07_RepeatMasker.sh` (`05_parseRM.pl`), `08_plot_div.sh/.R` | Bin RepeatMasker divergence and plot the TE landscape (relative insertion age) | `plots/08_diversity/` |
| **2. Gene annotation** | `09_run_maker.sh` | MAKER with Trinity transcripts, TAIR10 + UniProt proteins, EDTA TE library, Augustus (*arabidopsis* model) | MAKER datastore* |
| | `10_merge.sh` | Merge per-contig GFF3/FASTA | `output/09_MAKER/assembly.all.maker.*` |
| | `11_rename.sh`, `12_consistent_ids.sh` | Copy and rename gene IDs (prefix `ALT`) | `output/09_MAKER/final/` |
| | `13_InterProScan.sh`, `14_update.sh` | Pfam domains + GO terms, added to the GFF3 | `*.renamed.iprscan.gff` |
| | `15_extrct_filt.sh` | AED distribution, Pfam-by-AED, quality filter (AED < 1 and/or Pfam), filtered GFF3/FASTA | `filtered.genes.renamed.gff3` |
| | `16_Busco.sh` | Longest isoform per gene, BUSCO (brassicales_odb10) on proteins + transcripts, AGAT statistics | `busco_*/`, `annotation.stat` |
| **3. Function** | `17_run_blastp.sh` | BLASTP vs. UniProt (Viridiplantae, reviewed) and TAIR10, best hits, UniProt functions added to GFF3/FASTA, summary numbers | `output/comparisons/blastp_summary.txt` |
| **4. Comparative genomics** | `18_prepare_genespace.sh` | BED + longest-protein FASTA for Mh-0, TAIR10 and three Lian et al. accessions, ID validation | `output/comparisons/GENESPACE/{bed,peptide}` |
| | `20_run_genespace.sh` (`19_genespace.R`) | GENESPACE (OrthoFinder + MCScanX) and pan-gene table | `GENESPACE/results`, `pangenes/` |
| | `21_run_process_pangenome.sh/.R` | Core / dispensable / accession-specific entries, copy-variable core entries (all vs. syntenic genes) | `GENESPACE/process_pangenome/` |
| | `22_run_riparian.sh/.R` | Riparian synteny plots (Mh-0 vs. TAIR10, all five genomes) | `GENESPACE/processed_riparian/` |
| **Report** | `23_report_figures.sh/.R` | `Additional plots for report (e.g. extended TE superfamilies)` | `report/graphics/` |

\* `09_run_maker.sh` writes to `output/MAKER`; this folder was renamed to `output/09_MAKER` before running step 10.

## Notes

- **Gene ID prefix:** IDs carry the course-template prefix `ALT` (e.g. `ALT0026100`), but all genes belong to Mh-0.
- **BLAST ID fix (step 17):** BLAST+ 2.15 reports bare UniProt accessions (`Q9LX89`), whereas `maker_functional_fasta/gff`
  look up the full header ID (`sp|Q9LX89|FB255_ARATH`). Without remapping, `maker_functional_fasta` silently labels every
  protein "Protein of unknown function" and `maker_functional_gff` crashes (`Can't use string ("") as a HASH ref`).
  Step 17 therefore maps the best hits back to full IDs before calling both tools.
- **RNA-seq evidence** (Trinity transcripts from the assembly course) does not originate from Mh-0, so splice sites
  are mainly supported by protein homology.
- `05_parseRM.pl` is taken from [4ureliek/Parsing-RepeatMaskerOutputs](https://github.com/4ureliek/Parsing-RepeatMaskerOutputs).