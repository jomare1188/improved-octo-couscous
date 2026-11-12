#!/usr/bin/env Rscript
# H2: does 8008 remember more than 5000? (allele-group level, tillering mother leaves)
#  1. Joint DESeq2 model of all leaves of both genotypes (~ 0 + genotype_condition), so the two genotypes
#     can be compared directly. Interaction contrasts = 8008 effect minus 5000 effect:
#       int_D1_vs_W  (1st drought response), int_D3_vs_W (3rd), int_D3_vs_D1 (memory; 5000 Mt_D1 is one plant).
#  2. Which genotype amplifies its drought response more at D3? Drought genes are chosen from the average response
#     of both genotypes and both cycles (symmetric, so the choice does not favour either genotype or cycle); then
#     cameraPR tests the interaction statistic on the induced and on the repressed set.
#  3. Design check: 5000 has 1 Mt_D1 plant and 2 Mt_W plants. 8008 is refitted with the same design
#     (every choice of 1 of its 3 D1 plants and 2 of its 3 W plants, 9 fits, same model as 03_de.R)
#     to see how many D3 vs D1 DEGs and memory genes a 1-plant D1 group produces by itself.
# Run from analysis/:  qsub -N h2 code/run_r.sh code/10_genotype_cycle.R
source("code/common.R")
suppressMessages({ library(DESeq2); library(limma); library(BiocParallel) })
register(SerialParam())
OUT <- file.path(AN, "results/h2_genotype_cycle"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta <- read_meta(); txi <- readRDS(file.path(AN, "data/txi_group.rds"))

# ---------------------------------------------------------------- 1. joint model
leaf <- meta[meta$set != "bud", ]
leaf$grp <- factor(paste0("g", leaf$genotype, "_", leaf$condition))
dds <- DESeqDataSetFromTximport(subset_txi(txi, leaf$sample), leaf, ~ 0 + grp)
dds <- DESeq(dds, quiet = TRUE)
saveRDS(dds, file.path(OUT, "dds_joint_leaf.rds"))
rn <- resultsNames(dds)
con <- function(w) { v <- setNames(numeric(length(rn)), rn); k <- paste0("grp", names(w))
                     stopifnot(all(k %in% rn)); v[k] <- w; v }
G <- function(g, c) sprintf("g%s_%s", g, c)
cons <- list()
for (g in c("8008", "5000")) {
  cons[[paste0(g, "_Mt_D1_vs_Mt_W")]]  <- setNames(c(1, -1), c(G(g, "Mt_D1"), G(g, "Mt_W")))
  cons[[paste0(g, "_Mt_D3_vs_Mt_W")]]  <- setNames(c(1, -1), c(G(g, "Mt_D3"), G(g, "Mt_W")))
  cons[[paste0(g, "_Mt_D3_vs_Mt_D1")]] <- setNames(c(1, -1), c(G(g, "Mt_D3"), G(g, "Mt_D1")))
}
for (k in c("Mt_D1_vs_Mt_W", "Mt_D3_vs_Mt_W", "Mt_D3_vs_Mt_D1")) {
  a <- cons[[paste0("8008_", k)]]; b <- cons[[paste0("5000_", k)]]
  cons[[paste0("int_", k)]] <- c(a, -b)
}
cons[["avg_drought"]] <- setNames(rep(c(0.25, 0.25, -0.5), 2),
  c(G("8008", "Mt_D1"), G("8008", "Mt_D3"), G("8008", "Mt_W"), G("5000", "Mt_D1"), G("5000", "Mt_D3"), G("5000", "Mt_W")))

status_of <- function(k) if (grepl("D1", k) && (grepl("^5000", k) || grepl("^int", k) || k == "avg_drought")) "exploratory (5000 Mt_D1 n=1)" else "conclusive"
res <- list(); sm <- list()
for (k in names(cons)) {
  r <- results(dds, contrast = con(cons[[k]]), alpha = PADJ)
  dt <- data.table(id = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                   stat = r$stat, pvalue = r$pvalue, padj = r$padj)
  dt[, deg := call_deg(padj, log2FC)]
  fwrite(dt, file.path(OUT, sprintf("joint_%s.tsv.gz", k)), sep = "\t")
  res[[k]] <- dt
  sm[[k]] <- data.table(contrast = k, status = status_of(k), tested = sum(!is.na(dt$padj)),
                        up = sum(dt$deg == "up"), down = sum(dt$deg == "down"),
                        padj05_any_lfc = sum(dt$padj < PADJ, na.rm = TRUE))
}
sm <- rbindlist(sm); fwrite(sm, file.path(OUT, "joint_deg_summary.tsv"), sep = "\t"); print(sm)

# ---------------------------------------------------------------- 2. who amplifies more at D3?
avg <- res[["avg_drought"]]
sets <- list(induced = avg[deg == "up", id], repressed = avg[deg == "down", id])
cat("drought sets (average response):", lengths(sets), "\n")
amp <- list(); cam <- list()
for (k in c("int_Mt_D3_vs_Mt_D1", "int_Mt_D1_vs_Mt_W", "int_Mt_D3_vs_Mt_W")) {
  x <- res[[k]][!is.na(stat)]
  idx <- lapply(sets, function(s) which(x$id %in% s))
  r <- cameraPR(setNames(x$stat, x$id), idx)
  cam[[k]] <- data.table(contrast = k, set = rownames(r), r)
}
cam <- rbindlist(cam)
# reading: induced set + Direction Up = the 8008 value is higher than the 5000 value
cam[, reading := fifelse((set == "induced") == (Direction == "Up"), "8008 stronger", "5000 stronger")]
fwrite(cam, file.path(OUT, "interaction_camera.tsv"), sep = "\t"); print(cam)

# per-genotype amplification on the same genes: sign(response) x (D3 - D1); > 0 = stronger at D3
for (s in names(sets)) for (g in c("8008", "5000")) {
  x <- res[[paste0(g, "_Mt_D3_vs_Mt_D1")]][id %in% sets[[s]]]
  sg <- if (s == "induced") 1 else -1
  amp[[paste(s, g)]] <- data.table(set = s, genotype = g, n = nrow(x),
    median_amplification = round(median(sg * x$log2FC, na.rm = TRUE), 3),
    share_stronger_at_D3 = round(mean(sg * x$log2FC > 0, na.rm = TRUE), 3),
    sig_stronger = sum(x$padj < PADJ & sg * x$log2FC > 0, na.rm = TRUE),
    sig_weaker = sum(x$padj < PADJ & sg * x$log2FC < 0, na.rm = TRUE))
}
amp <- rbindlist(amp)
ia <- res[["int_Mt_D3_vs_Mt_D1"]]
for (s in names(sets)) {
  sg <- if (s == "induced") 1 else -1; x <- ia[id %in% sets[[s]] & !is.na(padj)]
  amp <- rbind(amp, data.table(set = s, genotype = "interaction (8008 - 5000)", n = nrow(x),
    median_amplification = round(median(sg * x$log2FC), 3), share_stronger_at_D3 = round(mean(sg * x$log2FC > 0), 3),
    sig_stronger = sum(x$padj < PADJ & sg * x$log2FC > 0), sig_weaker = sum(x$padj < PADJ & sg * x$log2FC < 0)))
}
fwrite(amp, file.path(OUT, "amplification_by_genotype.tsv"), sep = "\t"); print(amp)
long_amp <- rbindlist(lapply(names(sets), function(s) rbindlist(lapply(c("8008", "5000"), function(g) {
  x <- res[[paste0(g, "_Mt_D3_vs_Mt_D1")]][id %in% sets[[s]]]
  data.table(set = s, genotype = g, id = x$id, amplification = (if (s == "induced") 1 else -1) * x$log2FC)
})))); fwrite(long_amp, file.path(OUT, "amplification_per_group.tsv.gz"), sep = "\t")

# ---------------------------------------------------------------- 3. 8008 refitted with the 5000 design
l8 <- meta[meta$genotype == "8008" & meta$set != "bud", ]
W8 <- l8$sample[l8$condition == "Mt_W"]; D18 <- l8$sample[l8$condition == "Mt_D1"]
fit_counts <- function(md, tag) {
  md$condition <- factor(md$condition)
  dd <- DESeq(DESeqDataSetFromTximport(subset_txi(txi, md$sample), md, ~ condition), quiet = TRUE)
  get <- function(a, b) { r <- results(dd, contrast = c("condition", a, b), alpha = PADJ)
    data.table(id = rownames(r), lfc = r$log2FoldChange, padj = r$padj)[, deg := call_deg(padj, lfc)] }
  a <- get("Mt_D1", "Mt_W"); b <- get("Mt_D3", "Mt_D1"); d3 <- get("Mt_D3", "Mt_W")
  sgn <- function(x) fifelse(x$deg == "up", "+", fifelse(x$deg == "down", "-", "="))
  cl <- paste0(sgn(a), sgn(b)); cl[sgn(a) == "="] <- "other"
  s <- a$deg != "" | d3$deg != ""
  data.table(design = tag, D1_vs_W = sum(a$deg != ""), D3_vs_W = sum(d3$deg != ""), D3_vs_D1 = sum(b$deg != ""),
             mem_pp = sum(cl == "++"), mem_mm = sum(cl == "--"), mem_pm = sum(cl == "+-"), mem_mp = sum(cl == "-+"),
             r_D1_D3 = round(cor(a$lfc[s], d3$lfc[s], use = "complete.obs"), 3),
             median_abs_D3_vs_D1_on_responsive = round(median(abs(b$lfc[s]), na.rm = TRUE), 3))
}
ds <- list(fit_counts(l8, "8008 full (W 3, D1 3, D3 2)"))
for (wd in W8) for (dk in D18) {
  keep <- setdiff(l8$sample, c(wd, setdiff(D18, dk)))
  tag <- sprintf("8008 as 5000 (W 2, D1 1 = %s)", sub(".*_(\\d+[a-z_]*_S\\d+)$", "\\1", dk))
  ds[[length(ds) + 1]] <- fit_counts(l8[keep, ], tag); cat(tag, "done\n")
}
ds <- rbindlist(ds)
# 5000 as fitted in 03_de.R, for reference
s5 <- fread(file.path(AN, "results/de/group/deg_summary.tsv"), colClasses = list(character = "genotype"))[genotype == "5000"]
mc5 <- fread(file.path(AN, "results/de/group/memory_classes_5000_tillering.tsv.gz"))
h1 <- fread(file.path(AN, "results/report_tables/h1_d1_vs_d3_concordance.tsv"), colClasses = list(character = "genotype"))[genotype == "5000"]
ds <- rbind(ds, data.table(design = "5000 actual (W 2, D1 1, D3 2)",
  D1_vs_W = s5[contrast == "Mt_D1_vs_Mt_W", up + down], D3_vs_W = s5[contrast == "Mt_D3_vs_Mt_W", up + down],
  D3_vs_D1 = s5[contrast == "Mt_D3_vs_Mt_D1", up + down],
  mem_pp = sum(mc5$class == "++"), mem_mm = sum(mc5$class == "--"), mem_pm = sum(mc5$class == "+-"), mem_mp = sum(mc5$class == "-+"),
  r_D1_D3 = h1$pearson_lfc,
  median_abs_D3_vs_D1_on_responsive = { m <- mc5[(!is.na(D1_vs_W_padj) & D1_vs_W_padj < PADJ & abs(D1_vs_W) >= LFC) |
                                              (!is.na(D3_vs_W_padj) & D3_vs_W_padj < PADJ & abs(D3_vs_W) >= LFC)]
                                         round(median(abs(m$D3_vs_D1), na.rm = TRUE), 3) }), fill = TRUE)
fwrite(ds, file.path(OUT, "design_check_8008_as_5000.tsv"), sep = "\t"); print(ds)

# ---------------------------------------------------------------- figures
P <- plot_setup()
ds[, who := fifelse(grepl("full", design), "8008, full design", fifelse(grepl("as 5000", design), "8008, 5000-like design", "5000"))]
lng <- melt(ds, id.vars = c("design", "who"), measure.vars = c("D3_vs_D1", "mem_pp", "mem_mm", "r_D1_D3"))
lng[, variable := factor(variable, c("D3_vs_D1", "mem_pp", "mem_mm", "r_D1_D3"),
                         c("DEGs, D3 vs D1", "memory genes ++", "memory genes --", "r, D1 vs D3 response"))]
col3 <- c("8008, full design" = P$geno[["8008"]], "8008, 5000-like design" = "#9DB9D3", "5000" = P$geno[["5000"]])
p1 <- ggplot(lng, aes(who, value, colour = who)) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 2.6) +
  facet_wrap(~ variable, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = col3, guide = "none") + labs(x = NULL, y = NULL,
  title = "Does a one-plant D1 group by itself inflate the D3 vs D1 changes?",
  subtitle = "8008 refitted with every choice of 1 of its 3 D1 plants and 2 of its 3 watered plants (the 5000 design), 9 fits") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p1, file.path(OUT, "design_check.png"), 11, 4.4)

la <- long_amp[!is.na(amplification)]
la[, set := factor(set, c("induced", "repressed"), c("drought-induced groups", "drought-repressed groups"))]
p2 <- ggplot(la, aes(genotype, pmax(pmin(amplification, 4), -4), fill = genotype)) +
  geom_hline(yintercept = 0, colour = P$muted, linewidth = 0.4) +
  geom_violin(colour = NA, alpha = 0.55, scale = "width") + geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
  facet_wrap(~ set) + scale_fill_manual(values = P$geno, guide = "none") +
  labs(x = "Genotype", y = "Amplification at D3 (log2FC, clamped to ±4)",
       title = "How much stronger is the 3rd drought response than the 1st?",
       subtitle = "Amplification = D3 vs D1 change in the direction of the drought response; > 0 = stronger at D3.\nGroups chosen from the average response of both genotypes and cycles (joint model).")
save_plot(p2, file.path(OUT, "amplification.png"), 8, 4.8)
cat("done\n")
