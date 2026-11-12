#!/usr/bin/env Rscript
# Figures (allele-group level):
#   1. PCA: all libraries, and one panel per genotype x sample set
#   2. UpSet of DEGs per genotype; every contrast is split into up and down sets
#   3. Top 20 GO BP terms per contrast, up and down in separate plots
# Colours come from scico, and were checked with code/palette_check.R (adjacent/all-pair OKLab dE under
# deutan/protan/tritan simulation, plus contrast on white):
#   treatment history  managua 0.25/0.45/0.80   blue = watered, orange = D1 / tillering, dark red = D3 / maturation
#   sample set         glasgow 0/0.35/0.70
#   up / down          vik 0.15 / 0.85          (diverging poles: brown-red = up, blue = down)
#   GO -log10(p)       lajolla (up), devon (down); a single ramp each, light = weaker, dark = stronger
suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
                   library(scico); library(SummarizedExperiment); library(matrixStats); library(stringr) })
D   <- "/Storage/data1/jorge.munoz/DOLORES/analysis"
OUT <- file.path(D, "results/plots"); dir.create(file.path(OUT, "go"), showWarnings = FALSE, recursive = TRUE)

ink <- "#222222"; muted <- "#6B6B6B"; grid_col <- "#E6E6E6"
theme_set(theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
        axis.text = element_text(colour = muted), axis.title = element_text(colour = ink),
        plot.title = element_text(colour = ink, face = "bold", size = 12),
        plot.subtitle = element_text(colour = muted, size = 9), plot.caption = element_text(colour = muted, size = 8),
        legend.title = element_text(colour = ink, size = 9), legend.text = element_text(colour = ink, size = 8),
        plot.background = element_rect(fill = "white", colour = NA)))

m3   <- scico(101, palette = "managua")[round(c(0.80, 0.25, 0.45) * 100) + 1]   # watered, D1/T, D3/M
hist_col <- c(watered = m3[1], first = m3[2], second = m3[3])
set_col  <- setNames(scico(101, palette = "glasgow")[round(c(0, 0.35, 0.70) * 100) + 1],
                     c("Bud", "Mother leaf", "Propagule leaf"))
ud <- scico(101, palette = "vik")[round(c(0.85, 0.15) * 100) + 1]; dir_col <- c(up = ud[1], down = ud[2], mixed = "#8C8C8C")

# ---------------------------------------------------------------- 1. PCA
meta <- fread(file.path(D, "data/samples.tsv"), colClasses = list(character = "genotype")); setDF(meta)
rownames(meta) <- meta$sample
excl <- fread(file.path(D, "data/excluded_samples.tsv"))$sample
vsd  <- readRDS(file.path(D, "data/vst_group.rds")); m <- assay(vsd)[, meta$sample]
meta$set_lab  <- c(bud = "Bud", mother_leaf = "Mother leaf", propagule_leaf = "Propagule leaf")[meta$set]
meta$short    <- sub("^[0-9]+_[A-Za-z]+_[A-Z0-9]+_", "", meta$sample)
meta$excluded <- meta$sample %in% excl
meta$watered  <- meta$condition %in% c("W", "Mt_W", "Mm_W", "Fw_W")
# treatment history -> 3 colour classes; second factor -> shape
meta$history <- with(meta, ifelse(set == "bud", c(W = "watered", T = "first", M = "second")[condition],
                         ifelse(set == "mother_leaf", ifelse(grepl("_W$", condition), "watered", ifelse(grepl("D1", condition), "first", "second")),
                                c(Fw = "watered", Ft = "first", Fm = "second")[substr(condition, 1, 2)])))

pca_df <- function(idx, ntop = 2000) {
  x <- m[, idx, drop = FALSE]; x <- x[order(rowVars(x), decreasing = TRUE)[seq_len(min(ntop, nrow(x)))], ]
  p <- prcomp(t(x)); ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  list(df = cbind(meta[idx, ], PC1 = p$x[, 1], PC2 = p$x[, 2]), ve = ve)
}
lab_layer <- function() geom_text_repel(aes(label = short), size = 2.2, colour = muted, max.overlaps = 40,
                                        segment.colour = grid_col, min.segment.length = 0.2, seed = 1)
excl_layer <- function(df) geom_point(data = df[df$excluded, ], shape = 4, size = 4.2, stroke = 0.9, colour = ink)

# all libraries: colour = sample set, shape = genotype, filled = drought-exposed / hollow = watered
a <- pca_df(seq_len(ncol(m)))
g_all <- ggplot(a$df, aes(PC1, PC2)) +
  geom_point(aes(colour = set_lab, shape = genotype, fill = ifelse(watered, "white", set_col[set_lab])), size = 3, stroke = 1) +
  scale_fill_identity() + scale_colour_manual(values = set_col, name = "Sample set") +
  scale_shape_manual(values = c("5000" = 21, "8008" = 24), name = "Genotype",
                     guide = guide_legend(override.aes = list(fill = muted, colour = muted))) +
  excl_layer(a$df) + lab_layer() +
  labs(title = "PCA of all 51 libraries", subtitle = "Filled = drought-exposed (or from stressed mothers), hollow = watered. x = excluded (genotype mismatch).",
       x = sprintf("PC1 (%.1f%%)", a$ve[1]), y = sprintf("PC2 (%.1f%%)", a$ve[2]),
       caption = "VST of allele-group counts, 2,000 most variable groups")
ggsave(file.path(OUT, "pca_all.png"), g_all, width = 9, height = 6.5, dpi = 300, bg = "white")

panel <- function(s, g) {
  # excluded libraries are left out of per-genotype panels: their alleles belong to the other genotype and drive the PCs
  i <- which(meta$set == s & meta$genotype == g & !meta$excluded); a <- pca_df(i); df <- a$df
  if (s == "bud") {
    df$hist_lab <- c(watered = "W", first = "T (tillering-stressed)", second = "M (maturation-stressed)")[df$history]
    df$shape_by <- "bud"; shp_scale <- scale_shape_manual(values = c(bud = 21), guide = "none")
  } else if (s == "mother_leaf") {
    df$hist_lab <- c(watered = "W", first = "D1", second = "D3")[df$history]
    df$shape_by <- ifelse(grepl("^Mt", df$condition), "Tillering", "Maturation")
    shp_scale <- scale_shape_manual(values = c(Tillering = 21, Maturation = 22), name = "Phase",
                                    guide = guide_legend(override.aes = list(fill = muted, colour = muted)))
  } else {
    df$hist_lab <- c(watered = "Fw (mother W)", first = "Ft (mother T)", second = "Fm (mother M)")[df$history]
    df$shape_by <- ifelse(df$watered, "Watered", "Drought")
    shp_scale <- scale_shape_manual(values = c(Watered = 1, Drought = 21), name = "Propagule",
                                    guide = guide_legend(override.aes = list(fill = muted, colour = muted)))
  }
  cols <- setNames(hist_col[df$history], df$hist_lab); cols <- cols[!duplicated(names(cols))]
  lev <- names(cols)[order(match(names(cols), c("W", "D1", "D3", "T (tillering-stressed)", "M (maturation-stressed)",
                                               "Fw (mother W)", "Ft (mother T)", "Fm (mother M)")))]
  df$hist_lab <- factor(df$hist_lab, levels = lev)
  ggplot(df, aes(PC1, PC2)) +
    geom_point(aes(shape = shape_by, fill = hist_lab, colour = hist_lab), size = 3, stroke = 1) +
    scale_fill_manual(values = cols, name = "Treatment") + scale_colour_manual(values = cols, name = "Treatment") + shp_scale +
    lab_layer() +
    labs(title = sprintf("%s · %s", g, c(bud = "Buds", mother_leaf = "Mother leaves", propagule_leaf = "Propagule leaves")[s]),
         x = sprintf("PC1 (%.1f%%)", a$ve[1]), y = sprintf("PC2 (%.1f%%)", a$ve[2]))
}
pp <- list()
for (g in c("5000", "8008")) for (s in c("bud", "mother_leaf", "propagule_leaf")) {
  p <- panel(s, g); pp[[paste(g, s)]] <- p
  ggsave(file.path(OUT, sprintf("pca_%s_%s.png", s, g)), p, width = 6.5, height = 5, dpi = 300, bg = "white")
}
grid <- wrap_plots(pp, ncol = 3, byrow = TRUE) +
  plot_annotation(title = "PCA per genotype and sample set", caption = "The 3 libraries with a genotype mismatch are left out of these panels. VST of allele-group counts, 2,000 most variable groups per panel.")
ggsave(file.path(OUT, "pca_panels.png"), grid, width = 19, height = 11, dpi = 300, bg = "white")
cat("PCA done\n")

# ---------------------------------------------------------------- 2. UpSet (built in ggplot)
sm <- fread(file.path(D, "results/de/group/deg_summary.tsv"), colClasses = list(character = "genotype"))
order_c <- c("T_vs_W", "M_vs_W", "T_vs_M", "Mt_D1_vs_Mt_W", "Mt_D3_vs_Mt_W", "Mt_D3_vs_Mt_D1",
             "Mm_D1_vs_Mm_W", "Mm_D3_vs_Mm_W", "Mm_D3_vs_Mm_D1", "Mm_W_vs_Mt_W",
             "Fw_D1_vs_Fw_W", "Ft_D1_vs_Fw_D1", "Fm_D1_vs_Fw_D1", "Ft_D1_vs_Fw_W", "Fm_D1_vs_Fw_W")
nice <- function(x) gsub("_vs_", " vs ", x)
tag  <- c(conclusive = "", exploratory = " [expl]", descriptive = " [descr]")

upset_plot <- function(g, n_top = 45) {
  s <- sm[genotype == g][order(match(contrast, order_c))]
  mem <- rbindlist(lapply(seq_len(nrow(s)), function(k) {
    dt <- fread(file.path(D, sprintf("results/de/group/de_%s_%s.tsv.gz", g, s$contrast[k])), select = c("id", "deg"))[deg != ""]
    dt[, set := paste0(nice(s$contrast[k]), tag[s$status[k]], ifelse(deg == "up", "  ▲ up", "  ▼ down"))]
    dt[, .(id, set, dir = deg)]
  }))
  set_lev <- unlist(lapply(seq_len(nrow(s)), function(k) paste0(nice(s$contrast[k]), tag[s$status[k]], c("  ▲ up", "  ▼ down"))))
  set_lev <- set_lev[set_lev %in% mem$set]
  key <- mem[, .(key = paste(sort(match(set, set_lev)), collapse = ","),
                 kdir = if (all(dir == "up")) "up" else if (all(dir == "down")) "down" else "mixed"), by = id]
  inter <- key[, .N, by = .(key, kdir)][order(-N)][seq_len(min(n_top, .N))]
  inter[, x := seq_len(.N)]
  shown <- sum(inter$N); total <- uniqueN(key$id)
  mat <- inter[, .(row = as.integer(strsplit(key, ",")[[1]])), by = .(x, key)]
  mat[, set := set_lev[row]][, dir := ifelse(grepl("▲", set), "up", "down")]
  allpts <- CJ(x = inter$x, set = set_lev)
  set_n <- mem[, .N, by = set][, set := factor(set, levels = rev(set_lev))][, dir := ifelse(grepl("▲", set), "up", "down")]
  yl <- factor(set_lev, levels = rev(set_lev))
  p_bar <- ggplot(inter, aes(x, N, fill = kdir)) + geom_col(width = 0.72) +
    geom_text(aes(label = N), vjust = -0.4, size = 2.1, colour = muted) +
    scale_fill_manual(values = dir_col, name = "Intersection", labels = c(up = "all up", down = "all down", mixed = "up and down")) +
    scale_x_continuous(expand = expansion(add = 0.6)) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = "Groups in intersection") +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "top")
  nset <- length(set_lev)
  # separators between contrasts (each contrast has up to 2 rows); rows are numbered from the bottom
  ctr_of <- sub("  [▲▼].*", "", rev(set_lev)); sep_y <- which(head(ctr_of, -1) != tail(ctr_of, -1)) + 0.5
  p_mat <- ggplot() +
    geom_hline(yintercept = sep_y, colour = grid_col, linewidth = 0.4) +
    geom_point(data = allpts, aes(x, factor(set, levels = levels(yl))), colour = "#E3E3E3", size = 1.8) +
    geom_line(data = mat, aes(x, factor(set, levels = levels(yl)), group = x), colour = muted, linewidth = 0.5) +
    geom_point(data = mat, aes(x, factor(set, levels = levels(yl)), colour = dir), size = 2.2) +
    scale_colour_manual(values = dir_col, guide = "none") +
    scale_y_discrete(limits = levels(yl)) +
    scale_x_continuous(expand = expansion(add = 0.6)) +
    labs(x = NULL, y = NULL) + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.text.y = element_text(colour = ink, size = 8))
  p_set <- ggplot(set_n, aes(N, set, fill = dir)) + geom_col(width = 0.7) +
    geom_hline(yintercept = sep_y, colour = grid_col, linewidth = 0.4) +
    scale_fill_manual(values = dir_col, guide = "none") + scale_x_reverse(expand = expansion(mult = c(0.12, 0))) +
    scale_y_discrete(limits = levels(yl)) +
    labs(x = "DEG groups", y = NULL) + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
  lay <- "#A\nCB"
  wrap_plots(A = p_bar, B = p_mat, C = p_set, design = lay, widths = c(1, 5), heights = c(1.1, 2)) +
    plot_annotation(title = sprintf("Shared DEGs across contrasts · genotype %s", g),
      subtitle = sprintf("%d largest intersections shown (%s of %s DEG groups). DEG = padj < 0.05, |log2FC| ≥ 1, allele-group level.",
                         nrow(inter), format(shown, big.mark = ","), format(total, big.mark = ",")),
      caption = "[expl] = a group has 1 replicate (informs, no conclusions). [descr] = maturation mother leaves (descriptive). Bud W vs stressed is confounded with sequencing run.")
}
for (g in c("5000", "8008")) {
  ggsave(file.path(OUT, sprintf("upset_degs_%s.png", g)), upset_plot(g), width = 20, height = 11, dpi = 300, bg = "white")
}
cat("UpSet done\n")

# ---------------------------------------------------------------- 3. GO BP top 20, up and down separately
go <- fread(file.path(D, "results/enrichment/go_enrichment.tsv"))[ontology == "BP"]
go[, c("comp", "dir") := tstrsplit(set, ":", fixed = TRUE)]
go <- go[!grepl("memory", comp)]
go[, fold := Significant / Expected][, logp := -log10(pmax(p, 1e-30))]
n_made <- 0
for (cp in unique(go$comp)) for (d in c("up", "down")) {
  x <- go[comp == cp & dir == d][order(p, -fold)][1:min(20, .N)]
  if (!nrow(x)) next
  g <- sub("_.*", "", cp); ctr <- sub("^[0-9]+_", "", cp); st <- sm[genotype == g & contrast == ctr]
  x[, label := str_wrap(Term, 48)][, label := factor(label, levels = rev(unique(label)))]
  # both ramps run light (weaker) -> dark (stronger); the light ends are trimmed so every dot keeps contrast on white
  pal <- if (d == "up") list(p = "lajolla", b = 0.15, e = 0.75, dir = -1) else list(p = "devon", b = 0.05, e = 0.62, dir = -1)
  ndeg <- if (d == "up") st$up else st$down
  pl <- ggplot(x, aes(fold, label)) +
    geom_segment(aes(x = 0, xend = fold, yend = label), colour = grid_col, linewidth = 0.6) +
    geom_point(aes(size = Significant, fill = logp), shape = 21, colour = "white", stroke = 0.8) +
    scale_fill_scico(palette = pal$p, begin = pal$b, end = pal$e, direction = pal$dir, name = expression(-log[10]~p)) +
    scale_size(range = c(2, 7.5), trans = "sqrt", name = "DEG groups\nin term",
               guide = guide_legend(override.aes = list(fill = muted))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = sprintf("%s · %s · %s-regulated", g, nice(ctr), d),
         subtitle = sprintf("Top %d GO biological process terms (topGO weight01, Fisher p < 0.01)\n%s DEG groups · %s",
                            nrow(x), format(ndeg, big.mark = ","), st$status),
         x = "Fold enrichment (observed / expected)", y = NULL) +
    theme(axis.text.y = element_text(colour = ink, size = 8.5, lineheight = 0.85), panel.grid.major.y = element_blank(),
          plot.title.position = "plot")
  ggsave(file.path(OUT, "go", sprintf("go_bp_%s_%s.png", cp, d)), pl, width = 8.5, height = 1.2 + 0.36 * nrow(x), dpi = 300, bg = "white", limitsize = FALSE)
  n_made <- n_made + 1
}
cat("GO plots:", n_made, "\n")
