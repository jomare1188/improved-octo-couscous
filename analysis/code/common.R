# Shared paths, sample sheet, DE readers and plot style for the follow-up scripts (10_ onwards).
# Source from analysis/:  source("code/common.R")
suppressMessages(library(data.table))
D      <- "/Storage/data1/jorge.munoz/DOLORES"
AN     <- file.path(D, "analysis")
NF     <- file.path(D, "nf-rnaseq/run_POJ2878")
ORTHO  <- file.path(D, "orthology/allele_groups")
PADJ <- 0.05; LFC <- 1

# Sample sheet without the excluded libraries, plus the sequencing run of each library.
read_meta <- function(keep_excluded = FALSE) {
  meta <- fread(file.path(AN, "data/samples.tsv"), colClasses = list(character = "genotype"))
  runs <- fread(file.path(AN, "results/qc/sequencing_runs.tsv"), colClasses = list(character = "run"))
  meta <- merge(meta, runs[, .(sample, run)], by = "sample", all.x = TRUE, sort = FALSE)
  if (!keep_excluded) meta <- meta[!sample %in% fread(file.path(AN, "data/excluded_samples.tsv"))$sample]
  setDF(meta); rownames(meta) <- meta$sample; meta
}

# DE tables written by 03_de.R
rd_de <- function(g, contrast, level = "group")
  fread(file.path(AN, "results/de", level, sprintf("de_%s_%s.tsv.gz", g, contrast)))
call_deg <- function(padj, lfc) fifelse(!is.na(padj) & padj < PADJ & abs(lfc) >= LFC, fifelse(lfc > 0, "up", "down"), "")

# tximport list restricted to some samples; features with >= 10 reads in >= min_n of them
subset_txi <- function(txi, idx, min_n = 2) {
  k <- rowSums(txi$counts[, idx, drop = FALSE] >= 10) >= min_n
  list(abundance = txi$abundance[k, idx], counts = txi$counts[k, idx], length = txi$length[k, idx],
       countsFromAbundance = "no")
}

# Group annotation: a readable label per allele group
group_labels <- function() {
  gi  <- fread(file.path(ORTHO, "group_info.tsv"))
  g2g <- fread(file.path(ORTHO, "gene2group.tsv"))
  ann <- fread(file.path(D, "annotation/go/gene_annotation.tsv"), select = c("gene_id", "preferred_name", "description"))
  a <- merge(g2g[, .(gene_id, group_id)], ann, by = "gene_id")
  top <- function(x) { x <- x[!is.na(x) & x != "" & x != "-"]; if (!length(x)) "" else names(sort(table(x), decreasing = TRUE))[1] }
  a <- a[, .(name = top(preferred_name), description = top(description)), by = group_id]
  merge(gi[, .(group_id, n_genes, Athaliana, Athaliana_symbols)], a, by = "group_id", all.x = TRUE)
}

# ---------------------------------------------------------------- plot style (same as 06_plots.R)
plot_setup <- function() {
  suppressMessages({ library(ggplot2); library(scico) })
  theme_set(theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.3),
          axis.text = element_text(colour = "#6B6B6B"), axis.title = element_text(colour = "#222222"),
          plot.title = element_text(colour = "#222222", face = "bold", size = 12),
          plot.subtitle = element_text(colour = "#6B6B6B", size = 9), plot.caption = element_text(colour = "#6B6B6B", size = 8),
          strip.text = element_text(colour = "#222222", face = "bold"),
          legend.title = element_text(colour = "#222222", size = 9), legend.text = element_text(colour = "#222222", size = 8),
          plot.background = element_rect(fill = "white", colour = NA)))
  m3 <- scico(101, palette = "managua")[round(c(0.80, 0.25, 0.45) * 100) + 1]
  ud <- scico(101, palette = "vik")[round(c(0.85, 0.15) * 100) + 1]
  # genotype pair: managua poles already validated for CVD (blue / orange)
  list(ink = "#222222", muted = "#6B6B6B", grid = "#E6E6E6",
       hist = c(watered = m3[1], first = m3[2], second = m3[3]),
       dir = c(up = ud[1], down = ud[2], mixed = "#8C8C8C"),
       geno = c("8008" = m3[1], "5000" = m3[2]),
       div = scico(101, palette = "vik"))
}
save_plot <- function(p, file, w, h) ggsave(file, p, width = w, height = h, dpi = 300, bg = "white")
