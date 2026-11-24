#!/usr/bin/env Rscript
# Publication figures for results/report_story.md, written to results/plots/story/ as PNG (600 dpi), SVG and PDF:
#   pca.*          PCA per cultivar x sample set (VST, 2,000 most variable groups per panel). The 3 libraries with
#                  the other cultivar's alleles are not used to fit the PCs; they are projected onto them (x).
#   upset_leaf.*   UpSet of the tillering mother-leaf DEGs (D1 vs W, D3 vs W, D3 vs D1), each contrast split
#                  into up and down sets; both cultivars share the same intersection columns.
# Palettes as in 06_plots.R (validated by code/palette_check.R). Run from analysis/ (code/run_r.sh or Rscript).
suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork); library(scico)
                   library(SummarizedExperiment); library(matrixStats) })
A   <- "/Storage/data1/jorge.munoz/DOLORES/analysis"
OUT <- file.path(A, "results/plots/story"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

fam <- "Noto Sans"
ink <- "#222222"; muted <- "#6B6B6B"; grid_col <- "#E6E6E6"; faint <- "#D9D9D9"
theme_set(theme_minimal(base_size = 8, base_family = fam) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = 0.25),
        axis.text = element_text(colour = muted, size = 6.5), axis.title = element_text(colour = ink, size = 7.5),
        plot.title = element_text(colour = ink, face = "bold", size = 8.5), plot.title.position = "plot",
        legend.title = element_text(colour = ink, size = 7, face = "bold"), legend.text = element_text(colour = ink, size = 6.5),
        legend.key.size = unit(3.2, "mm"), plot.tag = element_text(face = "bold", size = 10, colour = ink),
        plot.background = element_rect(fill = "white", colour = NA)))
m3 <- scico(101, palette = "managua")[round(c(0.80, 0.25, 0.45) * 100) + 1]          # watered, 1st/T/Ft, D3/M/Fm
hist_col <- c(watered = m3[1], first = m3[2], second = m3[3])
ud <- scico(101, palette = "vik")[round(c(0.85, 0.15) * 100) + 1]
dir_col <- c(up = ud[1], down = ud[2], mixed = "#8C8C8C")

save3 <- function(p, f, w, h) {                       # width/height in mm
  ggsave(file.path(OUT, paste0(f, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
  ggsave(file.path(OUT, paste0(f, ".svg")), p, width = w, height = h, units = "mm", bg = "white", device = svglite::svglite)
  ggsave(file.path(OUT, paste0(f, ".pdf")), p, width = w, height = h, units = "mm", bg = "white", device = cairo_pdf)
}

# ------------------------------------------------------------------------------------------- PCA
meta <- fread(file.path(A, "data/samples.tsv"), colClasses = list(character = "genotype")); setDF(meta)
rownames(meta) <- meta$sample
excl <- fread(file.path(A, "data/excluded_samples.tsv"))$sample
vsd  <- readRDS(file.path(A, "data/vst_group.rds")); m <- assay(vsd)[, meta$sample]
meta$excluded <- meta$sample %in% excl
meta$watered  <- meta$condition %in% c("W", "Mt_W", "Mm_W", "Fw_W")
meta$short    <- sub("^[0-9]+_[A-Za-z]+_[A-Z0-9]+_", "", meta$sample)
meta$history  <- with(meta, ifelse(set == "bud", c(W = "watered", T = "first", M = "second")[condition],
                         ifelse(set == "mother_leaf", ifelse(grepl("_W$", condition), "watered", ifelse(grepl("D1", condition), "first", "second")),
                                c(Fw = "watered", Ft = "first", Fm = "second")[substr(condition, 1, 2)])))
set_name <- c(bud = "Buds", mother_leaf = "Mother leaves", propagule_leaf = "Propagule leaves")
hist_lab <- list(bud = c(watered = "W (watered)", first = "T (stressed at tillering)", second = "M (stressed at maturation)"),
                 mother_leaf = c(watered = "W (watered)", first = "D1 (1st drought)", second = "D3 (3rd drought)"),
                 propagule_leaf = c(watered = "Fw (mother watered)", first = "Ft (mother T)", second = "Fm (mother M)"))
shape_def <- list(bud = NULL,
                  mother_leaf = list(name = "Phase", values = c(Tillering = 21, Maturation = 22)),
                  propagule_leaf = list(name = "Propagule", values = c(Watered = 1, Drought = 21)))

pca_panel <- function(s, g, ntop = 2000) {
  fit  <- which(meta$set == s & meta$genotype == g & !meta$excluded)
  proj <- which(meta$set == s & meta$genotype == g & meta$excluded)
  x <- m[, fit, drop = FALSE]; keep <- order(rowVars(x), decreasing = TRUE)[seq_len(ntop)]
  p <- prcomp(t(x[keep, ])); ve <- 100 * p$sdev^2 / sum(p$sdev^2)
  df <- cbind(meta[fit, ], PC1 = p$x[, 1], PC2 = p$x[, 2])
  if (length(proj)) {
    z <- scale(t(m[keep, proj, drop = FALSE]), center = p$center, scale = FALSE) %*% p$rotation[, 1:2]
    pj <- cbind(meta[proj, ], PC1 = z[, 1], PC2 = z[, 2])
  } else pj <- NULL
  df$hist <- factor(hist_lab[[s]][df$history], levels = hist_lab[[s]])
  df$shp  <- switch(s, bud = "Bud", mother_leaf = ifelse(grepl("^Mt", df$condition), "Tillering", "Maturation"),
                    propagule_leaf = ifelse(df$watered, "Watered", "Drought"))
  df$flag <- ifelse(df$sample %in% meta$sample[meta$short == "79_S6"] | df$short == "79_S6", "79_S6 (failed poly-A)", NA)
  cols <- setNames(hist_col, hist_lab[[s]])
  sd <- shape_def[[s]]
  shp_scale <- if (is.null(sd)) scale_shape_manual(values = c(Bud = 21), guide = "none") else
    scale_shape_manual(values = sd$values, name = sd$name, guide = guide_legend(order = 2, override.aes = list(fill = muted, colour = muted)))
  gp <- ggplot(df, aes(PC1, PC2)) +
    geom_point(aes(fill = hist, colour = hist, shape = shp), size = 2.1, stroke = 0.7) +
    geom_point(data = data.frame(hist = factor(levels(df$hist), levels = levels(df$hist)), PC1 = df$PC1[1], PC2 = df$PC2[1],
                                 shp = df$shp[1]), aes(fill = hist, colour = hist, shape = shp), alpha = 0, size = 2.1, stroke = 0.7, show.legend = TRUE) +  # keeps every key
    scale_fill_manual(values = cols, name = "Treatment", drop = FALSE, guide = guide_legend(order = 1, ncol = 1, override.aes = list(shape = 21, alpha = 1))) +
    scale_colour_manual(values = cols, name = "Treatment", drop = FALSE, guide = guide_legend(order = 1, ncol = 1, override.aes = list(shape = 21, alpha = 1))) +
    shp_scale +
    geom_text_repel(data = df[!is.na(df$flag), ], aes(label = flag), size = 2.1, colour = ink, family = fam,
                    segment.colour = muted, segment.size = 0.25, min.segment.length = 0, box.padding = 0.5, seed = 1) +
    scale_x_continuous(expand = expansion(mult = 0.12)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = sprintf("%s · %s", g, set_name[s]), x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]))
  if (!is.null(pj)) gp <- gp +
    geom_point(data = pj, aes(PC1, PC2), shape = 4, size = 2.2, stroke = 0.7, colour = ink, inherit.aes = FALSE) +
    geom_text_repel(data = pj, aes(PC1, PC2, label = paste0(short, " (excluded)")), inherit.aes = FALSE, size = 2.1,
                    colour = ink, family = fam, segment.colour = muted, segment.size = 0.25, min.segment.length = 0,
                    box.padding = 0.8, nudge_y = -0.12 * diff(range(df$PC2)), nudge_x = -0.12 * diff(range(df$PC1)), seed = 1)
  gp + theme(plot.margin = margin(4, 8, 4, 4))
}
cols_p <- lapply(names(set_name), function(s)
  (pca_panel(s, "8008") / pca_panel(s, "5000")) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "vertical", legend.box = "horizontal",
          legend.box.just = "top", legend.justification = "left"))
fig_pca <- wrap_plots(cols_p, ncol = 3) + plot_annotation(tag_levels = list(c("A", "D", "B", "E", "C", "F")))
save3(fig_pca, "pca", 183, 150)
cat("PCA done\n")

# ------------------------------------------------------------------------------------------- UpSet, tillering leaves
ctr <- c(Mt_D1_vs_Mt_W = "D1 vs W", Mt_D3_vs_Mt_W = "D3 vs W", Mt_D3_vs_Mt_D1 = "D3 vs D1")
set_lev <- as.vector(t(outer(ctr, c("up", "down"), paste)))            # D1 vs W up, D1 vs W down, ...
mem <- rbindlist(lapply(c("8008", "5000"), function(g) rbindlist(lapply(names(ctr), function(cn) {
  d <- fread(file.path(A, sprintf("results/de/group/de_%s_%s.tsv.gz", g, cn)), select = c("id", "deg"))[deg != ""]
  d[, .(g, id, set = paste(ctr[cn], deg), dir = deg)] }))))
key <- mem[, .(key = paste(sort(match(set, set_lev)), collapse = ","), deg_n = .N,
               kdir = if (all(dir == "up")) "up" else if (all(dir == "down")) "down" else "mixed"), by = .(g, id)]
inter <- key[, .N, by = .(g, key, kdir, deg_n)]
cols_x <- inter[, .(tot = sum(N), deg_n = deg_n[1], kdir = kdir[1]), by = key][order(deg_n, -tot)][, x := .I]
inter <- merge(CJ(g = c("8008", "5000"), key = cols_x$key), inter[, .(g, key, N)], all.x = TRUE)[is.na(N), N := 0L]
inter <- merge(inter, cols_x[, .(key, x, kdir)], by = "key")
fwrite(dcast(inter, x + key ~ g, value.var = "N")[, sets := sapply(strsplit(key, ","), function(k) paste(set_lev[as.integer(k)], collapse = " & "))][order(x), .(x, sets, `8008`, `5000`)],
       file.path(A, "results/report_tables/upset_leaf_intersections.tsv"), sep = "\t")
mat <- cols_x[, .(row = as.integer(strsplit(key, ",")[[1]])), by = .(x, key)][, set := set_lev[row]]
mat[, dir := ifelse(grepl("up$", set), "up", "down")]
yl <- rev(set_lev); nx <- nrow(cols_x)
deg_sep <- cols_x[, max(x), by = deg_n]$V1; deg_sep <- head(deg_sep, -1) + 0.5
sep_y <- c(2.5, 4.5)
blk <- cols_x[, .(mid = mean(range(x))), by = deg_n][, lab := c("in one set only", "in two sets", "in three sets")[deg_n]]

upset_one <- function(gg, note, legend = FALSE) {
  it <- inter[g == gg]
  setn <- mem[g == gg, .N, by = set][, dir := ifelse(grepl("up$", set), "up", "down")]
  ymax <- max(it$N) * 1.22
  p_bar <- ggplot(it, aes(x, N, fill = kdir)) +
    geom_vline(xintercept = deg_sep, colour = grid_col, linewidth = 0.3) +
    geom_col(width = 0.72) +
    geom_text(aes(label = ifelse(N > 0, format(N, big.mark = ","), "")), vjust = -0.45, size = 1.9, colour = muted, family = fam) +
    scale_fill_manual(values = dir_col, name = "Intersection", breaks = c("up", "down", "mixed"),
                      labels = c(up = "all up", down = "all down", mixed = "up in one, down in another")) +
    scale_x_continuous(limits = c(0.4, nx + 0.6), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0), labels = scales::label_comma()) +
    annotate("text", x = blk$mid, y = ymax, label = blk$lab, vjust = 1, size = 2.2, colour = muted, family = fam, fontface = "italic") +
    labs(x = NULL, y = "Groups", title = sprintf("%s%s", gg, note)) +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = if (legend) "left" else "none")
  p_mat <- ggplot() +
    geom_hline(yintercept = sep_y, colour = grid_col, linewidth = 0.3) +
    geom_point(data = CJ(x = seq_len(nx), set = set_lev), aes(x, set), colour = faint, size = 1.3) +
    geom_line(data = mat, aes(x, set, group = x), colour = muted, linewidth = 0.4) +
    geom_point(data = mat, aes(x, set, colour = dir), size = 1.6) +
    scale_colour_manual(values = dir_col, guide = "none") + scale_y_discrete(limits = yl) +
    scale_x_continuous(limits = c(0.4, nx + 0.6), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) + theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_blank())
  p_set <- ggplot(setn, aes(N, set, fill = dir)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = format(N, big.mark = ",")), hjust = 1.1, size = 1.9, colour = muted, family = fam) +
    scale_fill_manual(values = dir_col, guide = "none") +
    scale_x_reverse(limits = c(5600, 0), expand = c(0, 0), breaks = c(4000, 2000, 0), labels = scales::label_comma()) +
    scale_y_discrete(limits = yl, position = "right") +
    labs(x = "DEGs per set", y = NULL) +
    theme(panel.grid.major.y = element_blank(), axis.text.y.right = element_text(colour = ink, size = 6.8, hjust = 0.5,
          margin = margin(l = 2, r = 2)))
  if (legend) wrap_plots(A = p_bar, B = p_mat, C = p_set, D = guide_area(), design = "DA\nCB", widths = c(1.15, 4), heights = c(1.25, 1), guides = "collect")
  else wrap_plots(A = p_bar, B = p_mat, C = p_set, design = "#A\nCB", widths = c(1.15, 4), heights = c(1.25, 1))
}
fig_up <- (wrap_elements(full = upset_one("8008", "", TRUE)) / wrap_elements(full = upset_one("5000", "  (one D1 plant: the D1 contrasts are exploratory)"))) +
  plot_annotation(tag_levels = "A")
save3(fig_up, "upset_leaf", 180, 175)
cat("UpSet done\n")
