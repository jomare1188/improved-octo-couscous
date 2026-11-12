#!/usr/bin/env Rscript
# Where do the three genotype-mismatched libraries come from?
# Their genotype is clear from allele presence/absence (qc/genotype_check). Here: within the true genotype,
# which group does each look like?
#   FIDH35000_S29  labelled 5000 T bud, carries 8008 alleles -> 8008 T or 8008 M bud?
#   21DH8_S24      labelled 8008 M bud, carries 5000 alleles -> 5000 T or 5000 M bud?
#   (a tube swap between the two predicts 8008 M and 5000 T; a genotype-only mix-up predicts 8008 T and 5000 M)
#   251_S3         labelled 5000 Mt_W leaf, carries 8008 alleles -> which 8008 leaf group?
# Method: nearest centroid on VST allele-group expression. Genes = the 1,000 groups that best separate the
# candidate groups (ratio of between- to within-group variance), chosen without the library being scored.
# Score for group g = correlation of (library - mean of centroids) with (centroid g - mean of centroids).
# Calibration: the same rule, leave-one-out, on every correctly labelled library of the reference set.
# Run from analysis/:  qsub -N relabel code/run_r.sh code/19_relabel_check.R
source("code/common.R")
suppressMessages({ library(SummarizedExperiment); library(ggplot2) })
OUT <- file.path(AN, "results/qc/relabel_check"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
P <- plot_setup()
meta <- read_meta(keep_excluded = TRUE)
excl <- fread(file.path(AN, "data/excluded_samples.tsv"))$sample
vst <- assay(readRDS(file.path(AN, "data/vst_group.rds")))[, meta$sample]

flagged <- list(
  list(lib = "5000_Tillering_D_FIDH35000_S29", genotype = "8008", set = "bud", classes = c("T", "M")),
  list(lib = "8008_Maturation_D_21DH8_S24", genotype = "5000", set = "bud", classes = c("T", "M")),
  list(lib = "5000_Mt_W_251_S3", genotype = "8008", set = c("mother_leaf", "propagule_leaf"),
       classes = c("Mt_W", "Mt_D1", "Mt_D3", "Mm_W", "Fw_W", "Fw_D1", "Ft_D1")))

score <- function(x, ref, cls, ntop = 1000) {
  ref <- ref[rowMeans(ref) > 4, , drop = FALSE]                 # expressed groups
  ok <- names(which(table(cls) >= 1)); keep <- cls %in% ok; ref <- ref[, keep, drop = FALSE]; cls <- cls[keep]
  cen <- sapply(ok, function(k) rowMeans(ref[, cls == k, drop = FALSE]))
  within <- sapply(seq_len(nrow(ref)), function(i) {
    r <- ref[i, ] - cen[i, cls]; if (length(r) > length(ok)) sum(r^2) / (length(r) - length(ok)) else NA })
  between <- apply(cen, 1, var)
  f <- between / pmax(ifelse(is.na(within), median(within, na.rm = TRUE), within), 0.05)
  g <- names(sort(f, decreasing = TRUE))[seq_len(min(ntop, length(f)))]
  mu <- rowMeans(cen[g, , drop = FALSE])
  sapply(ok, function(k) cor(x[g] - mu, cen[g, k] - mu))
}

res <- list(); loo <- list(); pcs <- list()
for (fl in flagged) {
  ref_s <- meta$sample[meta$genotype == fl$genotype & meta$set %in% fl$set & meta$condition %in% fl$classes &
                       !meta$sample %in% excl]
  cls <- setNames(meta[ref_s, "condition"], ref_s)
  # leave-one-out on the correctly labelled libraries
  for (s in ref_s) {
    others <- setdiff(ref_s, s)
    sc <- score(vst[, s], vst[rownames(vst), others], cls[others])
    srt <- sort(sc, decreasing = TRUE)
    loo[[length(loo) + 1]] <- data.table(reference = paste(fl$genotype, paste(fl$set, collapse = "+")), sample = s,
      true = cls[s], predicted = names(srt)[1], correct = names(srt)[1] == cls[s],
      score_true = round(sc[cls[s]], 3), margin = round(srt[1] - srt[2], 3),
      true_class_still_present = cls[s] %in% cls[others])
  }
  sc <- score(vst[, fl$lib], vst[, ref_s], cls)
  srt <- sort(sc, decreasing = TRUE)
  res[[length(res) + 1]] <- data.table(sample = fl$lib, labelled = paste(meta[fl$lib, "genotype"], meta[fl$lib, "condition"]),
    true_genotype = fl$genotype, best = names(srt)[1], margin = round(srt[1] - srt[2], 3),
    scores = paste(sprintf("%s %.2f", names(srt), srt), collapse = "; "))
  # nearest libraries (correlation on the 2,000 most variable groups of the reference set)
  x <- vst[, c(ref_s, fl$lib)]; x <- x[order(apply(x, 1, var), decreasing = TRUE)[1:2000], ]
  cc <- cor(x)[fl$lib, ref_s]
  res[[length(res)]]$nearest <- paste(sprintf("%s (%s) %.3f", sub("^[0-9]+_", "", names(sort(cc, decreasing = TRUE))[1:3]),
                                              cls[names(sort(cc, decreasing = TRUE))[1:3]], sort(cc, decreasing = TRUE)[1:3]), collapse = "; ")
  pc <- prcomp(t(x - rowMeans(x)))
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2))
  pcs[[length(pcs) + 1]] <- data.table(panel = sprintf("%s → true genotype %s (PC1 %d%%, PC2 %d%%)", sub("^[0-9]+_.*_", "", fl$lib), fl$genotype, ve[1], ve[2]),
    sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
    group = c(cls[ref_s], "flagged library"))
}
res <- rbindlist(res); loo <- rbindlist(loo)
fwrite(res, file.path(OUT, "flagged_assignment.tsv"), sep = "\t")
fwrite(loo, file.path(OUT, "loo_calibration.tsv"), sep = "\t")
print(res); print(loo[, .(n = .N, correct = sum(correct), median_margin = median(margin)), by = .(reference, true_class_still_present)])

pd <- rbindlist(pcs)
g <- ggplot(pd, aes(PC1, PC2, colour = group, shape = group == "flagged library")) +
  geom_point(size = 3) + facet_wrap(~ panel, scales = "free", ncol = 3) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8), guide = "none") +
  scale_colour_manual(values = c(setNames(scico::scico(7, palette = "batlow", end = 0.9),
                                          c("Mt_W", "Mt_D1", "Mt_D3", "Mm_W", "Fw_W", "Fw_D1", "Ft_D1")),
                                 T = unname(P$hist[2]), M = unname(P$hist[3]), `flagged library` = "black"), name = NULL) +
  labs(title = "Where do the three mismatched libraries sit within their true genotype?",
       subtitle = "PCA on the 2,000 most variable allele groups of each reference set; star = flagged library")
save_plot(g, file.path(OUT, "relabel_pca.png"), 15, 5)
cat("done\n")
