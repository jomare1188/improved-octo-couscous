#!/usr/bin/env Rscript
# Sampling hour and the memory genes. Each sampling lasted from early morning to past midday and the hour of
# each plant was not recorded, so two groups of a contrast can differ in time of day. What the clock genes say:
#  1. Internal time of each library from six core clock genes. Between early morning and midday the morning
#     gene CCA1/LHY falls and the day/evening genes (PRR73, PRR37, PRR59/95, GI, TOC1) rise. Per library:
#       morning score     = log2 CCA1/LHY
#       day/evening score = mean log2 of the five day/evening genes
#       clock index       = day/evening score - morning score        (higher = later internal time)
#     all centred within tissue (leaves, buds). Expression = VST of allele-group counts.
#  2. Per contrast: difference of the two scores between the groups. A later sampling hour moves them in
#     opposite directions (morning down, day/evening up).
#  3. The clock genes and the THI1 groups in the main DE tables (03_de.R), contrast by contrast.
# A model with the clock index as covariate (~ group + index) was tried and is not informative: plants of one
# group share their internal time (99% of the index variance is between groups), so the index cannot be
# separated from the group, and every adjusted test loses all power.
# Output: results/clock/.   Rscript code/23_clock_check.R   (from analysis/; about 2 minutes)
source("code/common.R")
suppressMessages(library(DESeq2))
OUT <- file.path(AN, "results/clock"); dir.create(OUT, showWarnings = FALSE)
P <- plot_setup()

clock <- data.table(
  id   = c("SB:SORBI_3007G047400", "SB:SORBI_3001G411400", "SB:SORBI_3006G057866", "SB:SORBI_3002G275100",
           "SB:SORBI_3003G040900", "SB:SORBI_3004G216700"),
  gene = c("CCA1/LHY", "PRR73", "PRR37", "PRR59/95", "GI", "TOC1"),
  part = c("morning", rep("day/evening", 5)))
de_genes <- clock$gene[clock$part == "day/evening"]

meta_all <- read_meta()
meta_all$tissue2 <- ifelse(meta_all$set == "bud", "bud", "leaf")
meta_all$grp <- factor(paste(meta_all$genotype, meta_all$condition, sep = "_"))
txi <- readRDS(file.path(AN, "data/txi_group.rds"))

# ---- 1. scores per library
lib <- rbindlist(lapply(c("leaf", "bud"), function(ts) {
  m <- meta_all[meta_all$tissue2 == ts, ]; m$grp <- droplevels(m$grp)
  dds <- DESeqDataSetFromTximport(subset_txi(txi, m$sample), m, ~ grp)
  v <- assay(vst(estimateSizeFactors(dds), blind = FALSE))
  stopifnot(all(clock$id %in% rownames(v)))
  ce <- v[clock$id, m$sample]; rownames(ce) <- clock$gene
  ctr <- function(x) as.numeric(scale(x, scale = FALSE))
  data.table(tissue = ts, sample = m$sample, library = m$library, genotype = m$genotype, condition = m$condition, run = m$run,
             morning = ctr(ce["CCA1/LHY", ]), day_evening = ctr(colMeans(ce[de_genes, ])),
             clock_index = ctr(colMeans(ce[de_genes, ]) - ce["CCA1/LHY", ]), t(ce))
}))
lev <- c("W", "T", "M", "Mt_W", "Mt_D1", "Mt_D3", "Mm_W", "Mm_D1", "Mm_D3", "Fw_W", "Fw_D1", "Ft_D1", "Fm_D1")
lib[, condition := factor(condition, levels = lev)]; setorder(lib, tissue, genotype, condition)
num <- c("morning", "day_evening", "clock_index", clock$gene)
fwrite(copy(lib)[, (num) := lapply(.SD, round, 2), .SDcols = num], file.path(OUT, "clock_index.tsv"), sep = "\t")
gm <- lib[, .(n = .N, morning = mean(morning), day_evening = mean(day_evening), clock_index = mean(clock_index),
              index_min = min(clock_index), index_max = max(clock_index)), by = .(tissue, genotype, condition)]
fwrite(copy(gm)[, 5:9 := lapply(.SD, round, 2), .SDcols = 5:9], file.path(OUT, "clock_index_by_group.tsv"), sep = "\t")
for (ts in c("leaf", "bud")) {
  a <- anova(lm(clock_index ~ paste(genotype, condition), lib[tissue == ts]))
  cat(sprintf("%s: %.0f%% of the clock index variance is between groups; SD within groups %.2f; range of group means %.2f\n",
              ts, 100 * a$`Sum Sq`[1] / sum(a$`Sum Sq`), sqrt(a$`Mean Sq`[2]), diff(range(gm[tissue == ts, clock_index]))))
}

# ---- 2. per contrast
cons <- rbind(
  data.table(tissue = "leaf", a = c("Mt_D1", "Mt_D3", "Mt_D3", "Fw_D1", "Ft_D1", "Fm_D1"), b = c("Mt_W", "Mt_W", "Mt_D1", "Fw_W", "Fw_D1", "Fw_D1")),
  data.table(tissue = "bud",  a = c("T", "T", "M"), b = c("M", "W", "W")))
ct <- rbindlist(lapply(c("8008", "5000"), function(g) rbindlist(lapply(seq_len(nrow(cons)), function(i) {
  x <- gm[tissue == cons$tissue[i] & genotype == g]; A <- x[condition == cons$a[i]]; B <- x[condition == cons$b[i]]
  if (!nrow(A) || !nrow(B)) return(NULL)
  dm <- A$morning - B$morning; de <- A$day_evening - B$day_evening
  data.table(genotype = g, contrast = sprintf("%s_vs_%s", cons$a[i], cons$b[i]), n1 = A$n, n2 = B$n,
             d_morning = round(dm, 2), d_day_evening = round(de, 2), d_clock_index = round(de - dm, 2),
             pattern = fcase(abs(dm) < 0.5 & abs(de) < 0.5, "same internal time",
                             dm <= -0.5 & de >= 0.5, "first group later (both scores)",
                             dm >= 0.5 & de <= -0.5, "first group earlier (both scores)",
                             abs(dm) >= 0.5 & abs(de) < 0.5, "CCA1/LHY alone differs",
                             abs(dm) < 0.5 & abs(de) >= 0.5, "day/evening genes alone differ",
                             default = "both scores move the same way (not a time shift)"))
}))))
fwrite(ct, file.path(OUT, "clock_by_contrast.tsv"), sep = "\t"); print(ct)

# ---- 3. clock genes and THI1 in the main DE tables
lab <- group_labels()
thi <- lab[name == "THI1", .(id = group_id, gene = "THI1", part = "clock-controlled output")]
gl <- rbind(clock, thi)
dt <- rbindlist(lapply(c("8008", "5000"), function(g) rbindlist(lapply(c("Mt_D1_vs_Mt_W", "Mt_D3_vs_Mt_W", "Mt_D3_vs_Mt_D1"), function(cn) {
  d <- rd_de(g, cn)[id %in% gl$id, .(id, baseMean = round(baseMean), log2FC = round(log2FC, 2), padj = signif(padj, 2))]
  d[, `:=`(genotype = g, contrast = cn)]
}))))
dt <- merge(gl, dt, by = "id")
mc <- rbindlist(lapply(c("8008", "5000"), function(g)
  fread(file.path(AN, sprintf("results/de/group/memory_classes_%s_tillering.tsv.gz", g)))[id %in% gl$id, .(genotype = g, id, class)]))
w <- dcast(dt, genotype + gene + part + id + baseMean ~ contrast, value.var = c("log2FC", "padj"))
w <- merge(w, mc, by = c("genotype", "id"), all.x = TRUE)
setorder(w, -genotype, part, gene, -baseMean)
fwrite(w, file.path(OUT, "clock_genes_de.tsv"), sep = "\t")
print(w[genotype == "8008", .(gene, id, baseMean, class, D1_W = log2FC_Mt_D1_vs_Mt_W, D3_D1 = log2FC_Mt_D3_vs_Mt_D1, padj_D3_D1 = padj_Mt_D3_vs_Mt_D1)])

# ---- figure: the two scores per library
pl <- melt(lib, id.vars = c("tissue", "genotype", "condition"), measure.vars = c("morning", "day_evening"), variable.name = "score")
pl[, score := factor(c(morning = "Morning score (CCA1/LHY)", day_evening = "Day/evening score (PRR73, PRR37, PRR59/95, GI, TOC1)")[as.character(score)],
                     levels = c("Morning score (CCA1/LHY)", "Day/evening score (PRR73, PRR37, PRR59/95, GI, TOC1)"))]
pl[, kind := fcase(condition %in% c("W", "Mt_W", "Mm_W", "Fw_W"), "watered", condition %in% c("T", "Mt_D3", "Mm_D3"), "3rd drought / T",
                   default = "1st drought / M")]
pl[, panel := factor(paste(genotype, ifelse(tissue == "bud", "buds", "leaves")), levels = c("8008 leaves", "5000 leaves", "8008 buds", "5000 buds"))]
p <- ggplot(pl, aes(condition, value, fill = kind)) +
  geom_hline(yintercept = 0, colour = P$grid) +
  geom_point(shape = 21, colour = "white", size = 2.8, stroke = 0.5, position = position_jitter(width = 0.1, height = 0, seed = 1)) +
  scale_fill_manual(values = c("watered" = unname(P$hist[1]), "1st drought / M" = unname(P$hist[2]), "3rd drought / T" = unname(P$hist[3])), name = NULL) +
  facet_grid(score ~ panel, scales = "free_x", space = "free_x", labeller = label_wrap_gen(28)) +
  labs(title = "Internal time of each library, from its clock genes",
       subtitle = "A later sampling hour lowers the morning score and raises the day/evening score. Dots = plants; values are log2, centred within tissue.",
       x = NULL, y = "Score (log2, centred)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = P$ink, size = 8), legend.position = "bottom")
save_plot(p, file.path(OUT, "clock_scores.png"), 11, 6.2)
