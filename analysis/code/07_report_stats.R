#!/usr/bin/env Rscript
# Extra statistics quoted in results/report.md (allele-group level).
#   h1_d1_vs_d3_concordance.tsv  how similar the 1st and 3rd drought responses of tillering leaves are
#   h5_overlaps.tsv              overlap of DEG lists between tissues (Fisher test, universe = groups tested in both)
#   h4_propagule_response.tsv    size of the propagule drought response on genes that respond in mother leaves
# Run from analysis/: Rscript code/07_report_stats.R
suppressMessages(library(data.table))

DE  <- "results/de/group"
OUT <- "results/report_tables"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(g, c) fread(file.path(DE, sprintf("de_%s_%s.tsv.gz", g, c)))[, .(id, lfc = log2FC, padj, deg)]

# H1: D1 vs W against D3 vs W, on groups that are DEGs in either contrast
h1 <- rbindlist(lapply(c("8008", "5000"), function(g) {
  m <- merge(rd(g, "Mt_D1_vs_Mt_W"), rd(g, "Mt_D3_vs_Mt_W"), by = "id", suffixes = c(".d1", ".d3"))
  s <- m[deg.d1 != "" | deg.d3 != ""]
  both <- s[deg.d1 != "" & deg.d3 != ""]
  data.table(genotype = g, deg_D1 = sum(m$deg.d1 != ""), deg_D3 = sum(m$deg.d3 != ""), deg_either = nrow(s),
             shared_same_dir = sum(both$deg.d1 == both$deg.d3), shared_opposite = sum(both$deg.d1 != both$deg.d3),
             pearson_lfc = round(cor(s$lfc.d1, s$lfc.d3, use = "complete.obs"), 3),
             median_abs_lfc_D1 = round(median(abs(s$lfc.d1), na.rm = TRUE), 2),
             median_abs_lfc_D3 = round(median(abs(s$lfc.d3), na.rm = TRUE), 2))
}))
fwrite(h1, file.path(OUT, "h1_d1_vs_d3_concordance.tsv"), sep = "\t")

# H5: do the same groups change in mother leaves, buds and propagules?
ov <- function(g, c1, c2) {
  m <- merge(rd(g, c1)[!is.na(padj)], rd(g, c2)[!is.na(padj)], by = "id")
  x <- m$deg.x != ""; y <- m$deg.y != ""
  ft <- fisher.test(table(factor(x, c(FALSE, TRUE)), factor(y, c(FALSE, TRUE))))
  both <- m[x & y]
  data.table(genotype = g, A = c1, B = c2, universe = nrow(m), nA = sum(x), nB = sum(y), overlap = sum(x & y),
             odds_ratio = round(unname(ft$estimate), 2), p = signif(ft$p.value, 2),
             same_dir = sum(both$deg.x == both$deg.y), opposite_dir = sum(both$deg.x != both$deg.y))
}
pairs <- list(
  c("Mt_D1_vs_Mt_W", "T_vs_M"), c("Mt_D3_vs_Mt_D1", "T_vs_M"), c("Mt_D1_vs_Mt_W", "T_vs_W"),
  c("T_vs_W", "Ft_D1_vs_Fw_D1"), c("T_vs_M", "Ft_D1_vs_Fw_D1"), c("Mt_D3_vs_Mt_D1", "Ft_D1_vs_Fw_D1"),
  c("Mt_D1_vs_Mt_W", "Ft_D1_vs_Fw_D1"))
h5 <- rbindlist(c(
  lapply(pairs, function(p) rbind(ov("8008", p[1], p[2]), ov("5000", p[1], p[2]))),
  list(ov("5000", "M_vs_W", "Fm_D1_vs_Fw_D1"), ov("5000", "Mt_D1_vs_Mt_W", "Fm_D1_vs_Fw_D1"))))
h5[, fdr := signif(p.adjust(p, "BH"), 2)]
fwrite(h5, file.path(OUT, "h5_overlaps.tsv"), sep = "\t")

# H4: response size in propagules, on groups chosen independently from mother leaves (Mt_D1_vs_Mt_W DEGs).
# Choosing the genes from other plants avoids the bias of picking genes from the contrast being measured.
h4 <- rbindlist(lapply(c("8008", "5000"), function(g) {
  ref <- rd(g, "Mt_D1_vs_Mt_W")[deg != "", .(id, dir = deg, lfc_mother = lfc)]
  cs <- c("Fw_D1_vs_Fw_W", "Ft_D1_vs_Fw_W", if (g == "5000") "Fm_D1_vs_Fw_W")
  rbindlist(lapply(cs, function(c) {
    x <- merge(ref, rd(g, c), by = "id")
    data.table(genotype = g, contrast = c,
               n_mother_up = sum(x$dir == "up"), median_lfc_on_mother_up = round(median(x[dir == "up", lfc], na.rm = TRUE), 2),
               n_mother_down = sum(x$dir == "down"), median_lfc_on_mother_down = round(median(x[dir == "down", lfc], na.rm = TRUE), 2),
               slope_vs_mother = round(unname(coef(lm(lfc ~ 0 + lfc_mother, data = x))[1]), 2))
  }))
}))
fwrite(h4, file.path(OUT, "h4_propagule_response.tsv"), sep = "\t")

cat("written to", OUT, "\n")
