#!/usr/bin/env Rscript
# Explanatory figures for results/report_story.md (allele-group level), written to results/plots/story/:
#   design_runs.png        libraries per group and sequencing run (why some contrasts are confounded)
#   memory_classes.png     schematic of the memory classes (Ding et al.)
#   d1_vs_d3.png           1st vs 3rd drought response of tillering leaves, coloured by memory class
#   example_genes.png      per-plant expression of four example genes in tillering leaves
#   propagule_response.png how strongly propagules respond to drought, by mother history
#   bud_pal.png            PAL (phenylalanine ammonia-lyase) family expression in buds
# Palettes are the ones validated for 06_plots.R (code/palette_check.R): managua for treatment history,
# vik poles for up/down. Run from analysis/ (or via code/run_story.sh).
suppressMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork); library(scico) })
D    <- "/Storage/data1/jorge.munoz/DOLORES"
A    <- file.path(D, "analysis")
DE   <- file.path(A, "results/de/group")
OUT  <- file.path(A, "results/plots/story"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ink <- "#222222"; muted <- "#6B6B6B"; grid_col <- "#E6E6E6"; faint <- "#CFCFCF"
theme_set(theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
        axis.text = element_text(colour = muted), axis.title = element_text(colour = ink),
        plot.title = element_text(colour = ink, face = "bold", size = 12), plot.title.position = "plot",
        plot.subtitle = element_text(colour = muted, size = 9), plot.caption = element_text(colour = muted, size = 8),
        strip.text = element_text(colour = ink, face = "bold", size = 9.5),
        legend.title = element_text(colour = ink, size = 9), legend.text = element_text(colour = ink, size = 8.5),
        plot.background = element_rect(fill = "white", colour = NA)))
m3 <- scico(101, palette = "managua")[round(c(0.80, 0.25, 0.45) * 100) + 1]   # watered, 1st / T, 2nd / M
ud <- scico(101, palette = "vik")[round(c(0.85, 0.15) * 100) + 1]             # up, down
up_col <- ud[1]; down_col <- ud[2]; rev_col <- "#4D4D4D"
save <- function(p, f, w, h) ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 300, bg = "white")

meta <- fread(file.path(A, "data/samples.tsv"), colClasses = list(character = "genotype"))
excl <- fread(file.path(A, "data/excluded_samples.tsv"))$sample
meta[, excluded := sample %in% excl]

# ------------------------------------------------------------ fig1: libraries per group and sequencing run
fi <- fread(file.path(D, "nf-rnaseq/run_POJ2878/integrity_check/fileinfo.tsv"))
# runs field: "@VH00942:44:AACV7GJM5:1=15268173;..." -> keep the run with most reads
fi[, run := sapply(strsplit(runs, ";"), function(r) {
  r <- r[nzchar(r)]; n <- as.numeric(sub(".*=", "", r)); sub("^@[^:]+:([0-9]+):.*", "\\1", r[which.max(n)]) })]
rt <- merge(meta[, .(sample, library, genotype, set, condition, excluded)], fi[, .(library = target, run)], by = "library", all.x = TRUE)
fwrite(rt[, .(sample, library, genotype, set, condition, run, excluded)], file.path(A, "results/qc/sequencing_runs.tsv"), sep = "\t")
lev_c <- c("W", "T", "M", "Mt_W", "Mt_D1", "Mt_D3", "Mm_W", "Mm_D1", "Mm_D3", "Fw_W", "Fw_D1", "Ft_D1", "Fm_D1")
rt[, set_lab := factor(c(bud = "Buds", mother_leaf = "Mother leaves", propagule_leaf = "Propagule leaves")[set],
                       levels = c("Buds", "Mother leaves", "Propagule leaves"))]
rt[, row := factor(paste(genotype, condition), levels = rev(as.vector(outer(c("5000", "8008"), lev_c, paste))))]
rt[, run_lab := factor(paste("run", run), levels = paste("run", sort(unique(as.integer(run)))))]
cnt <- rt[!(excluded), .N, by = .(set_lab, row, run_lab)]
ex  <- rt[(excluded), .N, by = .(set_lab, row, run_lab)]
f1 <- ggplot(cnt, aes(run_lab, row)) +
  geom_tile(fill = "#4A4A4A", colour = "white", linewidth = 1.2) +
  geom_text(aes(label = N), colour = "white", size = 3.3, fontface = "bold") +
  geom_point(data = ex, shape = 4, size = 3, stroke = 0.9, colour = ink, position = position_nudge(x = 0.33)) +
  facet_grid(set_lab ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Which group was sequenced on which run",
       subtitle = "Number = libraries used in the models. x = a library excluded for carrying the other genotype's alleles.\nWatered buds (W) sit alone on run 44 and every stressed bud on run 45, so 'stressed vs watered bud' is also 'run 45 vs run 44'.",
       x = "Sequencing run (instrument VH00942, from the FASTQ headers)", y = NULL) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(colour = ink, size = 9), axis.text.x = element_text(colour = ink),
        strip.text.y = element_text(angle = 0, hjust = 0))
save(f1, "design_runs.png", 9, 8)
cat("fig1 done\n")

# ------------------------------------------------------------ fig2: memory classes, schematic
sch <- rbindlist(list(
  data.table(cls = "++", D1 =  2, D3 =  4, txt = "memory: stronger induction at D3"),
  data.table(cls = "+-", D1 =  2, D3 =  0, txt = "revised: induced, then back down"),
  data.table(cls = "+=", D1 =  2, D3 =  2, txt = "no memory: same induction"),
  data.table(cls = "0+", D1 =  0, D3 =  2, txt = "late responder: only at D3"),
  data.table(cls = "--", D1 = -2, D3 = -4, txt = "memory: stronger repression at D3"),
  data.table(cls = "-+", D1 = -2, D3 =  0, txt = "revised: repressed, then back up"),
  data.table(cls = "-=", D1 = -2, D3 = -2, txt = "no memory: same repression"),
  data.table(cls = "0-", D1 =  0, D3 = -2, txt = "late responder: only at D3")))
sch[, panel := factor(sprintf("%s  %s", cls, txt), levels = sprintf("%s  %s", cls, txt))]
sch[, first := ifelse(D1 > 0, "up", ifelse(D1 < 0, "down", ifelse(D3 > 0, "up", "down")))]
sl <- melt(sch, id.vars = c("cls", "panel", "first", "txt"), measure.vars = c("D1", "D3"), variable.name = "step", value.name = "y")
sl <- rbind(sch[, .(cls, panel, first, txt, step = "W", y = 0)], sl)
sl[, step := factor(step, levels = c("W", "D1", "D3"))]
f2 <- ggplot(sl, aes(step, y, group = panel, colour = first)) +
  geom_hline(yintercept = 0, colour = faint, linewidth = 0.4) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  scale_colour_manual(values = c(up = up_col, down = down_col), guide = "none") +
  facet_wrap(~ panel, nrow = 2, labeller = label_wrap_gen(24)) +
  scale_y_continuous(breaks = NULL) +
  labs(title = "Memory classes: how a gene can behave over repeated droughts",
       subtitle = "First sign = D1 vs watered (W). Second sign = D3 vs D1. Idealised expression, not data.",
       x = NULL, y = "Expression") +
  theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(colour = ink))
save(f2, "memory_classes.png", 10, 5.2)
cat("fig2 done\n")

# ------------------------------------------------------------ fig3: D1 vs D3 response
gi <- fread(file.path(D, "orthology/allele_groups/group_info.tsv"))
lab_ids <- c("SB:SORBI_3009G116700" = "RAB18 dehydrin", "SB:SORBI_3003G149100" = "LEA3/SAG21",
             "SB:SORBI_3001G144900" = "seed-storage protein", "SB:SORBI_3002G384400" = "THI1",
             "OG:OG0019289" = "embryo-specific protein 3")
clamp <- 12
sc <- rbindlist(lapply(c("8008", "5000"), function(g) {
  m <- fread(file.path(DE, sprintf("memory_classes_%s_tillering.tsv.gz", g)))
  m <- m[(D1_vs_W_padj < 0.05 & abs(D1_vs_W) >= 1) | (D3_vs_W_padj < 0.05 & abs(D3_vs_W) >= 1)]
  r <- cor(m$D1_vs_W, m$D3_vs_W, use = "complete.obs")
  m[, genotype := sprintf("%s  (r = %.2f, %s groups)", g, r, format(.N, big.mark = ","))]
  m[, g := g]
}))
sc[, genotype := factor(genotype, levels = c(unique(genotype[g == "8008"]), unique(genotype[g == "5000"])))]
sc[, kind := fcase(class == "++", "++  stronger induction at D3", class == "--", "--  stronger repression at D3",
                   class %in% c("+-", "-+"), "+- / -+  reversed at D3", default = "same response or other")]
sc[, kind := factor(kind, levels = c("same response or other", "++  stronger induction at D3", "--  stronger repression at D3", "+- / -+  reversed at D3"))]
sc[, `:=`(x = pmax(pmin(D1_vs_W, clamp), -clamp), y = pmax(pmin(D3_vs_W, clamp), -clamp))]
sc[, label := ifelse(g == "8008" & id %in% names(lab_ids), lab_ids[id], NA_character_)]
f3 <- ggplot(sc, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, colour = muted, linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = grid_col) + geom_vline(xintercept = 0, colour = grid_col) +
  geom_point(data = sc[kind == "same response or other"], colour = faint, size = 0.7) +
  geom_point(data = sc[kind != "same response or other"], aes(colour = kind, shape = kind), size = 2) +
  geom_text_repel(aes(label = label), size = 3, colour = ink, min.segment.length = 0, seed = 2, box.padding = 0.6,
                  segment.colour = muted, na.rm = TRUE) +
  scale_colour_manual(values = c(up_col, down_col, rev_col), name = "Memory class", drop = FALSE,
                      limits = levels(sc$kind)[-1]) +
  scale_shape_manual(values = c(16, 16, 17), name = "Memory class", limits = levels(sc$kind)[-1]) +
  facet_wrap(~ genotype) + coord_equal() +
  labs(title = "Is the 3rd drought response different from the 1st?",
       subtitle = "Each dot is an allele group that is a DEG at D1 or at D3 (vs watered). On the dashed line the two responses are equal.\nGrey dots = same response (or no significant D3 vs D1 change). Values beyond ±12 are drawn at the edge.",
       x = "log2 fold change, 1st drought (D1 vs W)", y = "log2 fold change, 3rd drought (D3 vs W)",
       caption = "5000 has one D1 plant, so its panel is exploratory.") +
  theme(legend.position = "bottom")
save(f3, "d1_vs_d3.png", 11, 6.6)
cat("fig3 done\n")

# ------------------------------------------------------------ fig4: example genes, per plant
ex_ids <- c("SB:SORBI_3009G116700" = "RAB18 dehydrin",
            "SB:SORBI_3003G149100" = "LEA3/SAG21",
            "SB:SORBI_3002G172000|SORBI_3002G260300|SORBI_3002G360701|SORBI_3006G030500" = "HAI1 (PP2C)",
            "SB:SORBI_3002G384400" = "THI1")
eg <- rbindlist(lapply(c("8008", "5000"), function(g) {
  nc <- fread(file.path(DE, sprintf("normcounts_%s_leaf.tsv.gz", g)))[id %in% names(ex_ids)]
  l <- melt(nc, id.vars = "id", variable.name = "sample", value.name = "count")
  l <- merge(l, meta[, .(sample, condition)], by = "sample")[condition %in% c("Mt_W", "Mt_D1", "Mt_D3")]
  mc <- fread(file.path(DE, sprintf("memory_classes_%s_tillering.tsv.gz", g)))[id %in% names(ex_ids), .(id, class)]
  l <- merge(l, mc, by = "id"); l[, genotype := g]
}))
eg[, gene := factor(ex_ids[id], levels = ex_ids)]
eg[, step := factor(c(Mt_W = "W", Mt_D1 = "D1", Mt_D3 = "D3")[condition], levels = c("W", "D1", "D3"))]
eg[, y := log2(count + 1)][, genotype := factor(genotype, levels = c("8008", "5000"))]
mn <-eg[, .(y = mean(y), n = .N), by = .(gene, genotype, step, class)]
cl <- unique(eg[, .(gene, genotype, class)])
ymax <- eg[, .(top = max(y)), by = gene]
cl <- merge(cl, ymax, by = "gene")
f4 <- ggplot(eg, aes(step, y)) +
  geom_line(data = mn, aes(group = 1), colour = faint, linewidth = 0.8) +
  geom_errorbar(data = mn, aes(ymin = y, ymax = y, colour = step), width = 0.45, linewidth = 1.1) +
  geom_point(aes(fill = step), shape = 21, colour = "white", size = 3, stroke = 0.6,
             position = position_jitter(width = 0.08, height = 0, seed = 3)) +
  geom_text(data = cl, aes(x = 0.6, y = top * 1.08, label = paste("class", class)), hjust = 0, size = 3, colour = muted) +
  scale_fill_manual(values = setNames(m3, c("W", "D1", "D3")), guide = "none") +
  scale_colour_manual(values = setNames(m3, c("W", "D1", "D3")), guide = "none") +
  facet_grid(gene ~ genotype, scales = "free_y") +
  labs(title = "Four genes, plant by plant (tillering mother leaves)",
       subtitle = "Dots = plants, bar = group mean. In 8008, RAB18 and LEA3 keep rising at D3 ('++'), HAI1 jumps at D1 and stays ('+='),\nand THI1 drops further ('--'). In 5000 (a single D1 plant) the same genes fall in other classes.",
       x = NULL, y = "log2 (normalised count + 1)") +
  theme(axis.text.x = element_text(colour = ink), strip.text.y = element_text(angle = 0, hjust = 0))
save(f4, "example_genes.png", 8.5, 9)
cat("fig4 done\n")

# ------------------------------------------------------------ fig5: propagule response size
pr <- rbindlist(lapply(c("8008", "5000"), function(g) {
  ref <- fread(file.path(DE, sprintf("de_%s_Mt_D1_vs_Mt_W.tsv.gz", g)))[deg != "", .(id, dir = deg)]
  cs <- c(Fw = "Fw_D1_vs_Fw_W", Ft = "Ft_D1_vs_Fw_W", if (g == "5000") c(Fm = "Fm_D1_vs_Fw_W"))
  rbindlist(lapply(names(cs), function(k) {
    x <- merge(ref, fread(file.path(DE, sprintf("de_%s_%s.tsv.gz", g, cs[[k]])))[, .(id, lfc = log2FC)], by = "id")
    x[, `:=`(prop = k, genotype = g)]
  }))
}))
pr[, prop_lab := factor(c(Fw = "Fw\nwatered mother", Ft = "Ft\ntillering-stressed\nmother", Fm = "Fm\nmaturation-stressed\nmother")[prop],
                        levels = c("Fw\nwatered mother", "Ft\ntillering-stressed\nmother", "Fm\nmaturation-stressed\nmother"))]
pr[, dir_lab := factor(ifelse(dir == "up", "Genes drought INDUCES in mother leaves", "Genes drought REPRESSES in mother leaves"),
                       levels = c("Genes drought INDUCES in mother leaves", "Genes drought REPRESSES in mother leaves"))]
pr[, genotype := factor(genotype, levels = c("8008", "5000"))]
med <- pr[, .(med = median(lfc, na.rm = TRUE)), by = .(genotype, dir_lab, prop_lab)]
hist_cols <- setNames(m3, levels(pr$prop_lab))
f5 <- ggplot(pr, aes(prop_lab, lfc)) +
  geom_hline(yintercept = 0, colour = muted, linewidth = 0.4) +
  geom_boxplot(aes(fill = prop_lab), outlier.shape = NA, width = 0.55, colour = ink, linewidth = 0.35, alpha = 0.9) +
  geom_label(data = med, aes(y = med, label = sprintf("%+.2f", med)), size = 2.8, colour = ink, fill = "white",
             label.size = 0, nudge_x = 0.42) +
  scale_fill_manual(values = hist_cols, guide = "none") +
  coord_cartesian(ylim = c(-6, 6)) +
  facet_grid(genotype ~ dir_lab, scales = "free_x", space = "free_x") +
  labs(title = "Propagules from stressed mothers repress drought-repressed genes less",
       subtitle = "Propagule drought response (log2FC vs watered propagules of watered mothers, Fw_W) on genes chosen from the\nmother-leaf D1 response. Box = middle half of the genes, label = median. Every box shares the same Fw_W reference.",
       x = NULL, y = "log2 fold change in propagules (D1 vs Fw_W)",
       caption = "Exploratory: Fw_D1 is one plant in each genotype. 8008 has no Fm propagules.") +
  theme(axis.text.x = element_text(colour = ink, size = 8.5), panel.grid.major.x = element_blank(),
        strip.text.y = element_text(angle = 0))
save(f5, "propagule_response.png", 10, 7.5)
cat("fig5 done\n")

# ------------------------------------------------------------ fig6: PAL family in buds
g2g <- fread(file.path(D, "orthology/allele_groups/gene2group.tsv"), select = c("gene_id", "group_id"))
ann <- fread(file.path(D, "annotation/go/gene_annotation.tsv"), select = c("gene_id", "preferred_name", "description"))
pal_genes <- ann[description == "Phenylalanine ammonia-lyase" | grepl("^PAL[0-9]*$", preferred_name), gene_id]
pal_groups <- unique(g2g[gene_id %in% pal_genes, group_id])
pb <- rbindlist(lapply(c("8008", "5000"), function(g) {
  nc <- fread(file.path(DE, sprintf("normcounts_%s_bud.tsv.gz", g)))[id %in% pal_groups]
  s <- colSums(as.matrix(nc[, -1]))
  data.table(sample = names(s), total = s, genotype = g)
}))
pb <- merge(pb, meta[, .(sample, condition)], by = "sample")
pb[, step := factor(c(W = "W\nwatered", T = "T\ntillering-stressed\n(135 days ago)", M = "M\nmaturation-stressed\n(6 days ago)")[condition],
                    levels = c("W\nwatered", "T\ntillering-stressed\n(135 days ago)", "M\nmaturation-stressed\n(6 days ago)"))]
pb[, y := log2(total + 1)][, genotype := factor(genotype, levels = c("8008", "5000"))]
pm <- pb[, .(y = mean(y)), by = .(genotype, step)]
bud_cols <- setNames(m3, levels(pb$step))
f6 <- ggplot(pb, aes(step, y)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#F2F2F2") +
  annotate("text", x = 1, y = Inf, label = "run 44", vjust = 1.6, size = 2.8, colour = muted) +
  annotate("text", x = 2.5, y = Inf, label = "run 45", vjust = 1.6, size = 2.8, colour = muted) +
  geom_errorbar(data = pm, aes(ymin = y, ymax = y, colour = step), width = 0.45, linewidth = 1.1) +
  geom_point(aes(fill = step), shape = 21, colour = "white", size = 3.4, stroke = 0.6,
             position = position_jitter(width = 0.08, height = 0, seed = 4)) +
  scale_x_discrete() + scale_fill_manual(values = bud_cols, guide = "none") + scale_colour_manual(values = bud_cols, guide = "none") +
  facet_wrap(~ genotype) +
  labs(title = "PAL genes in buds: tillering vs maturation stress goes opposite ways in the two genotypes",
       subtitle = sprintf("Sum of normalised counts over the %d allele groups annotated as phenylalanine ammonia-lyase. Dots = plants, bar = mean.\nT and M share a run, so T vs M is a clean comparison; W was sequenced on another run (grey band).", length(pal_groups)),
       x = NULL, y = "log2 (summed normalised count + 1)") +
  theme(axis.text.x = element_text(colour = ink, size = 8.5), panel.grid.major.x = element_blank())
save(f6, "bud_pal.png", 9, 5.2)
cat("fig6 done\n")
