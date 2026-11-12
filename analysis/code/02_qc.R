#!/usr/bin/env Rscript
# Sample QC on the imported counts: library sizes, VST PCA (all samples, per set, per set x genotype)
# and a sample-sample correlation heatmap.   Rscript 02_qc.R <gene|group>
suppressMessages({ library(DESeq2); library(data.table); library(ggplot2); library(ggrepel)
                   library(pheatmap); library(matrixStats) })
level <- commandArgs(TRUE)[1]; if (is.na(level)) level <- "gene"
D   <- "/Storage/data1/jorge.munoz/DOLORES/analysis"
OUT <- file.path(D, "results/qc", level); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta <- fread(file.path(D, "data/samples.tsv"), colClasses = list(character = "genotype")); setDF(meta); rownames(meta) <- meta$sample
txi  <- readRDS(file.path(D, sprintf("data/txi_%s.rds", level)))
meta <- meta[colnames(txi$counts), ]

# library size and detected features
lib <- data.frame(sample = meta$sample, set = meta$set, genotype = meta$genotype, condition = meta$condition,
                  assigned_reads = round(colSums(txi$counts)),
                  features_ge10 = colSums(txi$counts >= 10))
fwrite(lib, file.path(OUT, "library_sizes.tsv"), sep = "\t")

keep <- rowSums(txi$counts >= 10) >= 2
txi_f <- lapply(txi[c("abundance", "counts", "length")], function(m) m[keep, ]); txi_f$countsFromAbundance <- "no"
dds <- DESeqDataSetFromTximport(txi_f, meta, ~ 1)
vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, file.path(D, sprintf("data/vst_%s.rds", level)))
m <- assay(vsd)

pca_plot <- function(idx, title, file, ntop = 2000) {
  x <- m[, idx, drop = FALSE]
  x <- x[order(rowVars(x), decreasing = TRUE)[seq_len(min(ntop, nrow(x)))], ]
  p <- prcomp(t(x)); ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  df <- cbind(meta[idx, ], PC1 = p$x[, 1], PC2 = p$x[, 2], PC3 = p$x[, 3])
  fwrite(df, sub(".png$", ".tsv", file), sep = "\t")
  g <- ggplot(df, aes(PC1, PC2, colour = condition, shape = genotype)) + geom_point(size = 3) +
    geom_text_repel(aes(label = sub("^[0-9]+_[A-Za-z]+_[A-Z0-9]+_", "", sample)), size = 2.3, show.legend = FALSE,
                    max.overlaps = 30) +
    labs(title = title, x = sprintf("PC1 (%s%%)", ve[1]), y = sprintf("PC2 (%s%%)", ve[2])) + theme_bw()
  if (length(unique(df$set)) > 1) g <- g + facet_null() + aes(colour = set)
  ggsave(file, g, width = 8, height = 6, dpi = 150)
}
pca_plot(seq_len(ncol(m)), "All 51 libraries", file.path(OUT, "pca_all.png"))
for (s in unique(meta$set)) {
  pca_plot(which(meta$set == s), s, file.path(OUT, sprintf("pca_%s.png", s)))
  for (g in unique(meta$genotype)) {
    i <- which(meta$set == s & meta$genotype == g)
    if (length(i) >= 4) pca_plot(i, paste(s, g), file.path(OUT, sprintf("pca_%s_%s.png", s, g)))
  }
}
leaf <- which(meta$set != "bud")
for (g in unique(meta$genotype)) {
  i <- intersect(leaf, which(meta$genotype == g))
  pca_plot(i, paste("all leaves", g), file.path(OUT, sprintf("pca_leaves_%s.png", g)))
}

top <- m[order(rowVars(m), decreasing = TRUE)[1:5000], ]
cc <- cor(top, method = "spearman")
ann <- meta[, c("set", "genotype", "condition")]
pheatmap(cc, annotation_col = ann, fontsize = 6, filename = file.path(OUT, "sample_correlation.png"),
         width = 14, height = 12)
cat("QC done:", level, nrow(m), "features after filter\n")
