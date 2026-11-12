#!/usr/bin/env Rscript
# H5 with rank-based tests (no DEG cut-off). Is it the same genes from mother leaf to bud to propagule?
#  1. Spearman correlation of the Wald statistic between every pair of contrasts, per genotype.
#  2. For chosen pairs A -> B: the 300 groups most induced and the 300 most repressed in A, tested as sets in the
#     ranking of B with cameraPR (rank-based, allows for correlation between genes); both directions.
#  3. Rank-rank overlap maps (RRHO, stratified as in RRHO2): for the top i groups of A and top j of B, from either
#     end of both rankings, -log10 p of the hypergeometric overlap.
#  New contrast: 5000 Ft_D1 vs Fm_D1 (propagules of tillering- vs maturation-stressed mothers, both under drought),
#  the propagule mirror of bud T vs M. 8008 has no Fm propagules.
# Run from analysis/:  qsub -N h5rank code/run_r.sh code/12_h5_rank.R
source("code/common.R")
suppressMessages({ library(DESeq2); library(limma) })
OUT <- file.path(AN, "results/h5_rank"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
P <- plot_setup()

dds5 <- readRDS(file.path(AN, "results/de/group/dds_5000_leaf.rds"))
r <- results(dds5, contrast = c("condition", "Ft_D1", "Fm_D1"), alpha = PADJ)
new5 <- data.table(id = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                   stat = r$stat, pvalue = r$pvalue, padj = r$padj)[, deg := call_deg(padj, log2FC)]
fwrite(new5, file.path(OUT, "de_5000_Ft_D1_vs_Fm_D1.tsv.gz"), sep = "\t")
cat("5000 Ft_D1 vs Fm_D1 (exploratory, Ft n=1):", sum(new5$deg == "up"), "up,", sum(new5$deg == "down"), "down\n")
rm(dds5); invisible(gc())

sm <- fread(file.path(AN, "results/de/group/deg_summary.tsv"), colClasses = list(character = "genotype"))
stats <- list()
for (g in c("8008", "5000")) {
  for (k in sm[genotype == g, contrast]) stats[[paste(g, k)]] <- rd_de(g, k)[!is.na(stat), .(id, stat)]
}
stats[["5000 Ft_D1_vs_Fm_D1"]] <- new5[!is.na(stat), .(id, stat)]
groups_of <- function(k) strsplit(k, "_vs_")[[1]]

# ---------------------------------------------------------------- 1. correlation matrices
order_c <- c("Mt_D1_vs_Mt_W", "Mt_D3_vs_Mt_W", "Mt_D3_vs_Mt_D1", "Mm_D1_vs_Mm_W", "Mm_D3_vs_Mm_W", "Mm_D3_vs_Mm_D1",
             "Mm_W_vs_Mt_W", "T_vs_W", "M_vs_W", "T_vs_M", "Fw_D1_vs_Fw_W", "Ft_D1_vs_Fw_W", "Fm_D1_vs_Fw_W",
             "Ft_D1_vs_Fw_D1", "Fm_D1_vs_Fw_D1", "Ft_D1_vs_Fm_D1")
cmat <- list()
for (g in c("8008", "5000")) {
  ks <- intersect(order_c, sub(paste0("^", g, " "), "", grep(paste0("^", g, " "), names(stats), value = TRUE)))
  wide <- Reduce(function(a, b) merge(a, b, by = "id", all = TRUE),
                 lapply(ks, function(k) setnames(copy(stats[[paste(g, k)]]), "stat", k)))
  cc <- cor(as.matrix(wide[, ..ks]), method = "spearman", use = "pairwise.complete.obs")
  long <- data.table(genotype = g, a = rep(ks, each = length(ks)), b = rep(ks, length(ks)), rho = as.vector(cc))
  long[, shared_group := mapply(function(x, y) length(intersect(groups_of(x), groups_of(y))) > 0, a, b)]
  cmat[[g]] <- long
}
cmat <- rbindlist(cmat); fwrite(cmat, file.path(OUT, "contrast_correlations.tsv"), sep = "\t")
pl <- copy(cmat)[, `:=`(a = factor(a, order_c), b = factor(b, rev(order_c)))]
pl[, lab := fifelse(as.integer(a) == length(order_c) + 1 - as.integer(b), "", sprintf("%.2f", rho))]
p1 <- ggplot(pl, aes(a, b, fill = rho)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab, colour = abs(rho) > 0.55), size = 2.1, show.legend = FALSE) +
  geom_point(data = pl[shared_group == TRUE & as.character(a) != as.character(b)], shape = 0, size = 5.2, colour = P$muted, stroke = 0.3) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = P$ink)) +
  scale_fill_gradientn(colours = rev(P$div), limits = c(-1, 1), name = "Spearman rho") +
  facet_wrap(~ genotype, scales = "free") +
  labs(x = NULL, y = NULL, title = "How similar are the contrasts? Rank correlation of the Wald statistic",
       subtitle = "Squares mark pairs that share a group of plants (their correlation is partly built in).") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 7), axis.text.y = element_text(size = 7), panel.grid = element_blank())
save_plot(p1, file.path(OUT, "contrast_correlations.png"), 15, 7.5)

# ---------------------------------------------------------------- 2. directional set tests between tissues
pairs <- rbindlist(list(
  data.table(g = "both", a = "Mt_D1_vs_Mt_W", b = "T_vs_M", step = "mother leaf -> bud"),
  data.table(g = "both", a = "Mt_D3_vs_Mt_D1", b = "T_vs_M", step = "mother leaf -> bud"),
  data.table(g = "both", a = "Mt_D1_vs_Mt_W", b = "T_vs_W", step = "mother leaf -> bud (run-confounded)"),
  data.table(g = "both", a = "Mt_D1_vs_Mt_W", b = "M_vs_W", step = "mother leaf -> bud (run-confounded)"),
  data.table(g = "both", a = "T_vs_M", b = "Ft_D1_vs_Fw_D1", step = "bud -> propagule"),
  data.table(g = "both", a = "T_vs_W", b = "Ft_D1_vs_Fw_D1", step = "bud -> propagule (run-confounded)"),
  data.table(g = "5000", a = "M_vs_W", b = "Fm_D1_vs_Fw_D1", step = "bud -> propagule (run-confounded)"),
  data.table(g = "5000", a = "T_vs_M", b = "Ft_D1_vs_Fm_D1", step = "bud -> propagule (matched)"),
  data.table(g = "both", a = "Mt_D1_vs_Mt_W", b = "Ft_D1_vs_Fw_D1", step = "mother leaf -> propagule"),
  data.table(g = "both", a = "Mt_D3_vs_Mt_D1", b = "Ft_D1_vs_Fw_D1", step = "mother leaf -> propagule"),
  data.table(g = "5000", a = "Mt_D1_vs_Mt_W", b = "Fm_D1_vs_Fw_D1", step = "mother leaf -> propagule"),
  data.table(g = "5000", a = "Mm_D1_vs_Mm_W", b = "Fm_D1_vs_Fw_D1", step = "mother leaf -> propagule")))
pairs <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  x <- pairs[i]; gs <- if (x$g == "both") c("8008", "5000") else x$g
  rbindlist(lapply(gs, function(g) data.table(ga = g, gb = g, a = x$a, b = x$b, step = x$step)))
}))
cross <- c("Mt_D1_vs_Mt_W", "Mt_D3_vs_Mt_W", "Mt_D3_vs_Mt_D1", "T_vs_M", "T_vs_W", "Ft_D1_vs_Fw_D1", "Fw_D1_vs_Fw_W")
pairs <- rbind(pairs, data.table(ga = "8008", gb = "5000", a = cross, b = cross, step = "same contrast, 8008 vs 5000"))
TOPN <- 300
settest <- function(sa, sb, from, to) {
  m <- merge(sa, sb, by = "id", suffixes = c(".a", ".b"))
  st <- setNames(m$stat.b, m$id)
  top <- list(top_induced = m[order(-stat.a)][1:TOPN, id], top_repressed = m[order(stat.a)][1:TOPN, id])
  cp <- cameraPR(st, lapply(top, function(s) which(names(st) %in% s)))
  data.table(from = from, to = to, set = rownames(cp), cp,
             median_stat_in_B = sapply(rownames(cp), function(s) round(median(st[top[[s]]]), 2)))
}
res <- list()
for (i in seq_len(nrow(pairs))) {
  x <- pairs[i]
  sa <- stats[[paste(x$ga, x$a)]]; sb <- stats[[paste(x$gb, x$b)]]
  if (is.null(sa) || is.null(sb)) next
  m <- merge(sa, sb, by = "id")
  rho <- cor(m$stat.x, m$stat.y, method = "spearman")
  fa <- paste(x$ga, x$a); fb <- paste(x$gb, x$b)
  res[[i]] <- rbind(settest(sa, sb, fa, fb), settest(sb, sa, fb, fa))[, `:=`(step = x$step, n_common = nrow(m), rho = round(rho, 3))]
}
res <- rbindlist(res); res[, FDR := p.adjust(PValue, "BH")]
setcolorder(res, c("step", "from", "to", "n_common", "rho", "set"))
fwrite(res, file.path(OUT, "h5_rank_tests.tsv"), sep = "\t"); print(res[, .(step, from, to, rho, set, Direction, PValue = signif(PValue, 2), FDR = signif(FDR, 2), median_stat_in_B)])

# ---------------------------------------------------------------- 3. RRHO maps
rrho <- function(a, b, step = 200) {
  m <- merge(a, b, by = "id"); N <- nrow(m); h <- floor(N / 2); nb <- floor(h / step); ni <- (1:nb) * step
  quad <- function(ra, rb, qa, qb) {
    ia <- ceiling(ra / step); ib <- ceiling(rb / step); k <- ia <= nb & ib <= nb
    tab <- matrix(tabulate(ia[k] + (ib[k] - 1) * nb, nb * nb), nb, nb)
    cum <- t(apply(apply(tab, 2, cumsum), 1, cumsum))
    I <- matrix(ni, nb, nb); J <- matrix(ni, nb, nb, byrow = TRUE)
    lp <- -phyper(cum - 1, I, N - I, J, lower.tail = FALSE, log.p = TRUE) / log(10)
    lp[cum < I * J / N] <- 0
    data.table(x = if (qa == "up") as.vector(I) / N else 1 - as.vector(I) / N,
               y = if (qb == "up") as.vector(J) / N else 1 - as.vector(J) / N,
               nlp = as.vector(lp), w = step / N, qa = qa, qb = qb, i = as.vector(I), j = as.vector(J), overlap = as.vector(cum))
  }
  rbind(quad(rank(-m$stat.x), rank(-m$stat.y), "up", "up"), quad(rank(m$stat.x), rank(m$stat.y), "down", "down"),
        quad(rank(-m$stat.x), rank(m$stat.y), "up", "down"), quad(rank(m$stat.x), rank(-m$stat.y), "down", "up"))
}
maps <- rbindlist(list(
  data.table(g = "8008", a = "Mt_D1_vs_Mt_W", b = "T_vs_M"), data.table(g = "5000", a = "Mt_D1_vs_Mt_W", b = "T_vs_M"),
  data.table(g = "8008", a = "Mt_D1_vs_Mt_W", b = "Ft_D1_vs_Fw_D1"), data.table(g = "5000", a = "Mt_D1_vs_Mt_W", b = "Ft_D1_vs_Fw_D1"),
  data.table(g = "8008", a = "T_vs_M", b = "Ft_D1_vs_Fw_D1"), data.table(g = "5000", a = "T_vs_M", b = "Ft_D1_vs_Fm_D1")))
rr <- rbindlist(lapply(seq_len(nrow(maps)), function(i) {
  mp <- maps[i]; d <- rrho(stats[[paste(mp$g, mp$a)]], stats[[paste(mp$g, mp$b)]])
  d[, panel := sprintf("%s: %s (x) vs %s (y)", mp$g, mp$a, mp$b)]
}))
rr_max <- rr[, .SD[which.max(nlp)], by = .(panel, qa, qb)][, .(panel, A = qa, B = qb, max_nlp = round(nlp, 1), top_i = i, top_j = j, overlap)]
fwrite(rr_max, file.path(OUT, "rrho_maxima.tsv"), sep = "\t"); print(rr_max)
p3 <- ggplot(rr, aes(x, y, fill = pmin(nlp, 50))) + geom_tile(aes(width = w, height = w)) +
  geom_vline(xintercept = 0.5, colour = "white", linewidth = 0.6) + geom_hline(yintercept = 0.5, colour = "white", linewidth = 0.6) +
  facet_wrap(~ panel, ncol = 2) +
  scale_fill_scico(palette = "lajolla", direction = -1, name = "-log10 p\n(capped at 50)") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) + scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
  coord_fixed() + labs(x = "first contrast (x): most up  \u2192  most down", y = "second contrast (y): most up  \u2192  most down", title = "Rank-rank overlap maps",
    subtitle = "Bottom-left = up in both, top-right = down in both, other corners = opposite directions.\nColour = how much more the top-ranked groups of the two lists overlap than chance.") +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 7), strip.text = element_text(size = 8))
save_plot(p3, file.path(OUT, "rrho_maps.png"), 9, 12)
cat("done\n")
