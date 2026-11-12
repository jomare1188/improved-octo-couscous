#!/usr/bin/env Rscript
# Allele-level look at the candidates. Sugarcane carries ~8-12 copies (alleles) of most genes; each allele group
# was analysed as one unit. Here the alleles inside a group are compared: does one copy carry the response
# (or the memory), or do all copies move together?
#  1. Per candidate group, the allele-level log2FC (results/de/gene) for D1 vs W, D3 vs D1 (tillering leaves) and
#     T vs M (buds), with the chromosome copy of each allele (ChrNNgK = chromosome NN, copy K).
#  2. Heterogeneity between alleles: Cochran's Q on the allele log2FCs (weights 1/SE^2) and I^2.
#     Genome-wide for every group that responds at D1 (group-level DEG) and has >= 2 tested alleles; then for
#     the memory groups on D3 vs D1.
#  3. Who carries the change: each allele's share of the group's change in normalised counts, compared with its
#     share of the group's expression ("all copies move together" = change share ~ expression share).
#  4. Can Salmon tell the alleles apart? From the 30 bootstrap replicates of each library: inferential relative
#     variance (InfRV) of each allele and of its group, and the correlation between alleles of one group across
#     bootstraps (negative = reads traded between near-identical copies).
#  5. Allele presence: which alleles are expressed in 8008, in 5000, or in both.
# Candidates: memory classes ++, --, +-, -+ of both genotypes, the named genes of the report, PAL groups.
# Run from analysis/:  qsub -N alleles code/run_r.sh code/13_alleles.R
source("code/common.R")
suppressMessages({ library(tximport) })
OUT <- file.path(AN, "results/alleles"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
P <- plot_setup(); meta <- read_meta()
MINBM <- 10

g2g <- fread(file.path(ORTHO, "gene2group.tsv"))[, .(gene_id, group_id)]
gff <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"{print $1\"\\t\"$9}' %s", file.path(NF, "results_pilot/genome/Sh_POJ2878_v1.gff3")),
             header = FALSE, sep = "\t", col.names = c("chr", "attr"))
gff[, gene_id := sub(";.*", "", sub("^ID=", "", attr))][, attr := NULL]
gff[, `:=`(hom = as.integer(sub("Chr(\\d+)g.*", "\\1", chr)), copy = as.integer(sub(".*g(\\d+)$", "\\1", chr)))]
g2g <- merge(g2g, gff, by = "gene_id", all.x = TRUE)
lab <- group_labels()

# ---------------------------------------------------------------- candidates
named <- data.table(group_id = c("SB:SORBI_3009G116700", "SB:SORBI_3003G149100", "SB:SORBI_3001G144900", "OG:OG0019289",
                                 "SB:SORBI_3002G384400", "SB:SORBI_3002G172000|SORBI_3002G260300|SORBI_3002G360701|SORBI_3006G030500"),
                    gene = c("RAB18 dehydrin", "LEA3/SAG21", "seed-storage protein", "embryo-specific protein 3", "THI1", "HAI1 (PP2C)"))
mc <- rbindlist(lapply(c("8008", "5000"), function(g) fread(file.path(AN, sprintf("results/de/group/memory_classes_%s_tillering.tsv.gz", g)))[, genotype := g]))
lhy <- lab[grepl("^(LHY|CCA1)$", name), group_id]
if (length(lhy)) named <- rbind(named, data.table(group_id = lhy, gene = paste0("CCA1/LHY", if (length(lhy) > 1) paste0(" ", seq_along(lhy)) else "")))
pal <- lab[grepl("^PAL", name) | grepl("Phenylalanine ammonia-lyase", description), group_id]
tvm <- rbindlist(lapply(c("8008", "5000"), function(g) rd_de(g, "T_vs_M")[id %in% pal & !is.na(stat)][order(-abs(stat))][1][, genotype := g]))
named <- rbind(named, data.table(group_id = unique(tvm$id), gene = paste0("PAL (top bud T vs M, ", tvm[!duplicated(id), genotype], ")")))
memory_groups <- mc[class %in% c("++", "--", "+-", "-+"), .(genotype, group_id = id, class)]
cand <- unique(c(named$group_id, memory_groups$group_id, pal))
cat("candidate groups:", length(cand), " named:", nrow(named), "\n")

# ---------------------------------------------------------------- allele-level tables
contr <- c(D1_vs_W = "Mt_D1_vs_Mt_W", D3_vs_D1 = "Mt_D3_vs_Mt_D1", D3_vs_W = "Mt_D3_vs_Mt_W", T_vs_M = "T_vs_M")
al <- list(); grp_lvl <- list()
for (g in c("8008", "5000")) for (k in names(contr)) {
  x <- rd_de(g, contr[[k]], "gene")[, .(gene_id = id, baseMean, log2FC, lfcSE, padj)]
  x <- merge(x, g2g, by = "gene_id"); x[, `:=`(genotype = g, contrast = k)]
  al[[paste(g, k)]] <- x
  grp_lvl[[paste(g, k)]] <- rd_de(g, contr[[k]])[, .(group_id = id, group_lfc = log2FC, group_padj = padj, group_deg = deg)][, `:=`(genotype = g, contrast = k)]
}
al <- rbindlist(al); grp_lvl <- rbindlist(grp_lvl)

# expression share and change share from allele-level normalised counts
share <- list()
for (g in c("8008", "5000")) for (tis in c("leaf", "bud")) {
  nc <- fread(file.path(AN, sprintf("results/de/gene/normcounts_%s_%s.tsv.gz", g, tis)))
  nc <- merge(g2g[group_id %in% cand, .(gene_id, group_id)], nc, by.x = "gene_id", by.y = "id")
  m <- meta[meta$genotype == g & meta$sample %in% names(nc), ]
  conds <- if (tis == "leaf") c("Mt_W", "Mt_D1", "Mt_D3") else c("T", "M", "W")
  for (cn in conds) nc[, (cn) := rowMeans(.SD), .SDcols = m$sample[m$condition == cn]]
  share[[paste(g, tis)]] <- nc[, c("gene_id", "group_id", conds), with = FALSE][, genotype := g]
}
sl <- merge(share[["8008 leaf"]], share[["8008 bud"]][, .(gene_id, T, M, W)], by = "gene_id", all.x = TRUE)
s5 <- merge(share[["5000 leaf"]], share[["5000 bud"]][, .(gene_id, T, M, W)], by = "gene_id", all.x = TRUE)
sh <- rbind(sl, s5)
sh[, expr_share_leaf := (Mt_W + Mt_D1 + Mt_D3) / sum(Mt_W + Mt_D1 + Mt_D3), by = .(genotype, group_id)]
sh[, change_share_D1_vs_W := (Mt_D1 - Mt_W) / sum(Mt_D1 - Mt_W), by = .(genotype, group_id)]
sh[, change_share_D3_vs_D1 := (Mt_D3 - Mt_D1) / sum(Mt_D3 - Mt_D1), by = .(genotype, group_id)]
sh[, expr_share_bud := (T + M) / sum(T + M, na.rm = TRUE), by = .(genotype, group_id)]
sh[, change_share_T_vs_M := (T - M) / sum(T - M, na.rm = TRUE), by = .(genotype, group_id)]

wide <- dcast(al[group_id %in% cand], genotype + group_id + gene_id + chr ~ contrast, value.var = c("baseMean", "log2FC", "lfcSE", "padj"))
cand_tab <- merge(wide, sh[, !c("Mt_W", "Mt_D1", "Mt_D3", "T", "M", "W")], by = c("genotype", "group_id", "gene_id"), all.x = TRUE)
cand_tab <- merge(cand_tab, lab[, .(group_id, name, description)], by = "group_id", all.x = TRUE)
fwrite(cand_tab[order(genotype, group_id, chr)], file.path(OUT, "candidate_alleles.tsv"), sep = "\t")

# ---------------------------------------------------------------- 2. heterogeneity between alleles
het <- function(x) {   # x: alleles of one group in one contrast
  x <- x[baseMean >= MINBM & is.finite(lfcSE) & lfcSE > 0 & !is.na(log2FC)]
  k <- nrow(x); if (k < 2) return(NULL)
  w <- 1 / x$lfcSE^2; b <- sum(w * x$log2FC) / sum(w); Q <- sum(w * (x$log2FC - b)^2)
  data.table(n_alleles = k, pooled_lfc = b, Q = Q, p_het = pchisq(Q, k - 1, lower.tail = FALSE), I2 = max(0, (Q - (k - 1)) / Q),
             n_sig_up = sum(x$padj < PADJ & x$log2FC > 0, na.rm = TRUE), n_sig_down = sum(x$padj < PADJ & x$log2FC < 0, na.rm = TRUE),
             lfc_range = diff(range(x$log2FC)))
}
H <- al[, het(.SD), by = .(genotype, contrast, group_id)]
H <- merge(H, grp_lvl, by = c("genotype", "contrast", "group_id"), all.x = TRUE)
H[, fdr_het := p.adjust(p_het, "BH"), by = .(genotype, contrast)]
fwrite(H[group_id %in% cand], file.path(OUT, "candidate_heterogeneity.tsv"), sep = "\t")
# genome-wide reference sets
resp <- H[contrast == "D1_vs_W" & group_deg != ""]
memo <- merge(H[contrast == "D3_vs_D1"], memory_groups, by = c("genotype", "group_id"))
nonmem <- merge(H[contrast == "D3_vs_D1"], mc[class %in% c("+=", "-="), .(genotype, group_id = id, class)], by = c("genotype", "group_id"))
gw <- rbind(resp[, .(set = "D1 response, all responsive groups", contrast, genotype, I2, fdr_het, n_sig_up, n_sig_down)],
            memo[, .(set = "D3 vs D1, memory groups (++ -- +- -+)", contrast, genotype, I2, fdr_het, n_sig_up, n_sig_down)],
            nonmem[, .(set = "D3 vs D1, non-memory groups (+= -=)", contrast, genotype, I2, fdr_het, n_sig_up, n_sig_down)])
gws <- gw[, .(groups = .N, median_I2 = round(median(I2), 2), pct_heterogeneous_fdr05 = round(100 * mean(fdr_het < 0.05), 1),
              pct_alleles_both_directions = round(100 * mean(n_sig_up > 0 & n_sig_down > 0), 1)), by = .(genotype, set)]
fwrite(gws, file.path(OUT, "heterogeneity_summary.tsv"), sep = "\t"); print(gws)

# ---------------------------------------------------------------- 3. who carries the change (memory groups)
carry <- merge(sh, memory_groups, by = c("genotype", "group_id"))
carry <- carry[, {
  o <- order(-change_share_D3_vs_D1); cs <- cumsum(change_share_D3_vs_D1[o])
  .(n_alleles = .N, top_allele_change_share = round(max(change_share_D3_vs_D1), 2),
    top_allele_expr_share = round(expr_share_leaf[o][1], 2),
    alleles_for_80pct_of_change = which(cs >= 0.8)[1],
    r_change_vs_expression = if (.N >= 3) round(cor(change_share_D3_vs_D1, expr_share_leaf), 2) else NA_real_)
}, by = .(genotype, group_id, class)]
carry <- merge(carry, lab[, .(group_id, name)], by = "group_id", all.x = TRUE)
fwrite(carry[order(genotype, -top_allele_change_share)], file.path(OUT, "memory_change_carriers.tsv"), sep = "\t")
print(carry[, .(groups = .N, median_alleles = median(n_alleles), median_top_change_share = median(top_allele_change_share, na.rm = TRUE),
                median_r_change_vs_expr = median(r_change_vs_expression, na.rm = TRUE)), by = genotype])

# ---------------------------------------------------------------- 4. bootstrap uncertainty (named groups + memory groups)
bt_groups <- unique(c(named$group_id, memory_groups$group_id))
t2g <- fread(file.path(NF, "results/salmon/salmon.merged.tx2gene.tsv"))[, .(transcript_id, gene_id)]
bt_genes <- g2g[group_id %in% bt_groups, gene_id]
bt_tx <- t2g[gene_id %in% bt_genes]
inf <- list(); bcor <- list()
for (s in meta$sample) {
  f <- file.path(NF, "results/salmon", s, "quant.sf")
  tx <- suppressMessages(tximport(f, type = "salmon", txOut = TRUE))
  reps <- tx$infReps[[1]]; rownames(reps) <- rownames(tx$counts)
  reps <- reps[intersect(bt_tx$transcript_id, rownames(reps)), , drop = FALSE]
  gmat <- rowsum(reps, bt_tx$gene_id[match(rownames(reps), bt_tx$transcript_id)])
  grpm <- rowsum(gmat, g2g$group_id[match(rownames(gmat), g2g$gene_id)])
  infrv <- function(m) { mu <- rowMeans(m); v <- matrixStats::rowVars(m); data.table(id = rownames(m), mean = mu, infrv = pmax(v - mu, 0) / pmax(mu, 1e-8) + 0.01) }
  inf[[s]] <- rbind(infrv(gmat)[, level := "allele"], infrv(grpm)[, level := "group"])[, sample := s]
  gg <- split(rownames(gmat), g2g$group_id[match(rownames(gmat), g2g$gene_id)])
  bcor[[s]] <- rbindlist(lapply(names(gg), function(G) {
    a <- gg[[G]]; a <- a[rowMeans(gmat[a, , drop = FALSE]) >= 20]; if (length(a) < 2) return(NULL)
    cm <- suppressWarnings(cor(t(gmat[a, , drop = FALSE]))); v <- cm[upper.tri(cm)]
    data.table(group_id = G, sample = s, pairs = length(v), min_cor = min(v, na.rm = TRUE), share_negative = mean(v < -0.3, na.rm = TRUE))
  }))
  cat(s, "\n")
}
inf <- rbindlist(inf); bcor <- rbindlist(bcor)
infs <- inf[mean >= 20, .(mean_infrv = mean(infrv), libraries = .N), by = .(level, id)]
fwrite(infs, file.path(OUT, "bootstrap_infrv.tsv"), sep = "\t")
fwrite(bcor[, .(min_cor = min(min_cor), median_share_negative = median(share_negative)), by = group_id],
       file.path(OUT, "bootstrap_allele_correlation.tsv"), sep = "\t")
print(infs[, .(features = .N, median_infrv = round(median(mean_infrv), 3), pct_infrv_gt_0.5 = round(100 * mean(mean_infrv > 0.5), 1)), by = level])

# ---------------------------------------------------------------- 5. allele presence by genotype (named groups)
pres <- dcast(sh[group_id %in% named$group_id, .(genotype, group_id, gene_id, on = (Mt_W + Mt_D1 + Mt_D3) / 3 >= 5)],
              group_id + gene_id ~ genotype, value.var = "on", fill = FALSE)
pres <- merge(pres, g2g[, .(gene_id, chr)], by = "gene_id")
pres[, expressed_in := fifelse(`8008` & `5000`, "both", fifelse(`8008`, "8008 only", fifelse(`5000`, "5000 only", "neither")))]
fwrite(merge(pres, named, by = "group_id"), file.path(OUT, "named_allele_presence.tsv"), sep = "\t")

# ---------------------------------------------------------------- figures
allele_lab <- function(x) x[, allele := sprintf("chr %s.%s  %s", hom, copy, sub("^POJ2878", "", gene_id))]
ccol <- c("1st drought vs watered" = P$hist[["first"]], "3rd vs 1st drought" = P$hist[["second"]], "bud T vs M" = P$hist[["watered"]])
forest <- function(d, file, title, w, h) {
  d <- allele_lab(copy(d)); d[, `:=`(lo = log2FC - 1.96 * lfcSE, hi = log2FC + 1.96 * lfcSE)]
  d[, genotype := factor(genotype, c("8008", "5000"))]
  p <- ggplot(d, aes(pmax(pmin(log2FC, 12), -12), allele, colour = contrast)) +
    geom_vline(xintercept = 0, colour = P$muted, linewidth = 0.3) +
    geom_linerange(aes(xmin = pmax(lo, -12), xmax = pmin(hi, 12)), alpha = 0.6) +
    geom_point(aes(size = pmax(share, 0.01))) +
    facet_grid(gene ~ genotype + contrast, scales = "free_y", space = "free_y", labeller = labeller(gene = label_wrap_gen(14))) +
    scale_colour_manual(values = ccol, guide = "none") + scale_size_area(max_size = 3.5, name = "share of the\ngroup's expression") +
    labs(x = "log2FC of the allele (95% interval; clamped to ±12)", y = NULL, title = title,
         subtitle = "Each row is one allele: chromosome.copy and gene number. Alleles with baseMean >= 10 in that genotype.") +
    theme(axis.text.y = element_text(size = 5.5), strip.text.y = element_text(angle = 0, size = 7),
          strip.text.x = element_text(size = 7), legend.position = "top")
  save_plot(p, file, w, h)
}
fp <- merge(al[group_id %in% named[!grepl("^PAL", gene), group_id] & contrast %in% c("D1_vs_W", "D3_vs_D1") & baseMean >= MINBM], named, by = "group_id")
fp <- merge(fp, sh[, .(genotype, gene_id, share = expr_share_leaf)], by = c("genotype", "gene_id"), all.x = TRUE)
fp[, contrast := factor(contrast, c("D1_vs_W", "D3_vs_D1"), c("1st drought vs watered", "3rd vs 1st drought"))]
forest(fp, file.path(OUT, "named_alleles.png"), "Do all copies of a gene respond together? Tillering mother leaves", 12, 13)
pp <- merge(al[group_id %in% tvm$id & contrast == "T_vs_M" & baseMean >= MINBM], named, by = "group_id")
pp <- merge(pp, sh[, .(genotype, gene_id, share = expr_share_bud)], by = c("genotype", "gene_id"), all.x = TRUE)
pp[, contrast := factor(contrast, "T_vs_M", "bud T vs M")]
forest(pp, file.path(OUT, "pal_alleles.png"), "PAL alleles in buds: tillering- vs maturation-stressed", 8, 6)

gw[, set := factor(set, unique(gw$set))]
p2 <- ggplot(gw, aes(set, I2, fill = set)) + geom_violin(colour = NA, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") + facet_wrap(~ genotype) +
  scale_fill_manual(values = unname(P$hist), guide = "none") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, 18)) +
  labs(x = NULL, y = expression(I^2 ~ "(share of allele differences beyond noise)"),
       title = "How different are the alleles of one group?",
       subtitle = "I2 = 0: the alleles agree within noise; I2 near 1: they respond differently.")
save_plot(p2, file.path(OUT, "allele_heterogeneity.png"), 9, 5)

ib <- copy(infs)
p3 <- ggplot(ib, aes(level, pmin(mean_infrv, 5), fill = level)) + geom_violin(colour = NA, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") + scale_y_log10() +
  scale_fill_manual(values = unname(P$hist[c("watered", "first")]), guide = "none") +
  labs(x = NULL, y = "Inferential relative variance (log scale)", title = "How sure is Salmon about each count?",
       subtitle = "From 30 bootstrap replicates per library; candidate and memory groups. Higher = less certain.")
save_plot(p3, file.path(OUT, "bootstrap_uncertainty.png"), 5.5, 4.8)
cat("done\n")
