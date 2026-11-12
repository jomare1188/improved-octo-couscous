#!/usr/bin/env Rscript
# Partial check of the bud run effect. Every watered bud (W) is on run 44 and every stressed bud (T, M) on run 45,
# so "stressed vs watered" cannot be separated from "run 45 vs run 44" by design. What can be checked:
#  1. Technical quality of each library (MultiQC: rRNA, mapping, duplication, GC, adapters, fragment length).
#  2. Run fingerprints: does a contrast's log2FC depend on group GC content or length? Contrasts that cross runs
#     are compared with contrasts inside one run set (T vs M buds; leaf D3 vs D1; propagule drought vs watered).
#  3. Genes that should be stable: Arabidopsis reference genes (Czechowski et al. 2005 and classic qPCR genes) and
#     empirical controls (groups that vary least across all leaf libraries of the genotype, which span five runs
#     and every treatment). If the run effect is small, they barely move in T vs W, as in T vs M.
#  4. RUVg (RUVSeq) with the empirical controls: estimate an unwanted factor and refit DESeq2 ~ W_1 + condition.
#     With complete confounding this cannot fully separate run from treatment; it shows how much of the
#     signal is explained by a factor that the control genes share.
#  5. Does the bud W-vs-stressed difference look like the leaf drought programme (genes chosen from mother leaves)?
# Run from analysis/:  qsub -N budbatch code/run_r.sh code/11_bud_batch.R
source("code/common.R")
suppressMessages({ library(DESeq2); library(limma); library(RUVSeq); library(Biostrings); library(BiocParallel) })
register(SerialParam())
OUT <- file.path(AN, "results/bud_batch"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta_all <- read_meta(keep_excluded = TRUE); meta <- read_meta()
P <- plot_setup()

# ---------------------------------------------------------------- 1. QC metrics per library
MQ <- file.path(NF, "results/multiqc/multiqc_report_data")
gs <- fread(file.path(MQ, "multiqc_general_stats.txt"))
gs[, key := sub("(_S[0-9]+)_[12]$", "\\1", sub(" Read [12]$", "", Sample))]
cols <- c(rRNA_pct = "bowtie2_rrna_removal-overall_alignment_rate", mapped_pct = "salmon-percent_mapped",
          duplication_pct = "fastqc_raw-percent_duplicates", GC_pct = "fastqc_filtered-percent_gc",
          adapter_pct = "trim_galore-tg_pct_with_adapter", quality_trimmed_pct = "trim_galore-tg_pct_quality_trimmed",
          raw_reads_M = "fastqc_raw-total_sequences", compatible_fragment_ratio = "salmon-compatible_fragment_ratio")
num_mean <- function(v) { v <- suppressWarnings(as.numeric(v)); if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE) }
qc <- gs[, lapply(.SD, num_mean), by = key, .SDcols = unname(cols)]; setnames(qc, unname(cols), names(cols))
sal <- fread(file.path(MQ, "multiqc_salmon.txt"))[, .(key = Sample, fragment_length = frag_length_mean)]
lib <- fread(file.path(AN, "results/qc/group/library_sizes.tsv"))[, .(key = sample, reads_on_genes_M = assigned_reads / 1e6,
                                                                       groups_detected = features_ge10)]
qc <- Reduce(function(a, b) merge(a, b, by = "key", all.x = TRUE), list(qc, sal, lib))
qc <- merge(data.table(meta_all)[, .(key = sample, genotype, set, condition, run)], qc, by = "key")
qc[, excluded := key %in% fread(file.path(AN, "data/excluded_samples.tsv"))$sample]
fwrite(qc, file.path(OUT, "qc_metrics.tsv"), sep = "\t")
mcols <- c(names(cols), "fragment_length", "reads_on_genes_M", "groups_detected")
bud <- qc[set == "bud" & !excluded]
qtest <- rbindlist(lapply(mcols, function(m) {
  a <- bud[run == "44", get(m)]; b <- bud[run == "45", get(m)]
  data.table(metric = m, median_run44_W = round(median(a), 2), median_run45_TM = round(median(b), 2),
             leaf_median = round(median(qc[set != "bud" & !excluded, get(m)], na.rm = TRUE), 2),
             wilcox_p = signif(suppressWarnings(wilcox.test(a, b)$p.value), 2))
}))
fwrite(qtest, file.path(OUT, "qc_run44_vs_run45.tsv"), sep = "\t"); print(qtest)
ql <- melt(qc[excluded == FALSE], id.vars = c("key", "genotype", "set", "condition", "run"), measure.vars = mcols)
ql[, grp := fifelse(set == "bud", paste0("bud ", condition), fifelse(set == "mother_leaf", "mother leaf", "propagule leaf"))]
ql[, run := factor(run, sort(unique(run)))]
gcol <- c("bud W" = P$hist[["watered"]], "bud T" = P$hist[["first"]], "bud M" = P$hist[["second"]],
          "mother leaf" = "#9C9C9C", "propagule leaf" = "#5E5E5E")
p1 <- ggplot(ql, aes(run, value, colour = grp, shape = genotype)) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 2), size = 2) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = gcol, name = NULL) + scale_shape_manual(values = c("5000" = 16, "8008" = 17), name = "Genotype") +
  labs(x = "Sequencing run", y = NULL, title = "Technical quality of every library, by sequencing run",
       subtitle = "Watered buds are all on run 44, stressed buds on run 45. Excluded libraries left out.")
save_plot(p1, file.path(OUT, "qc_metrics.png"), 12, 7)

# ---------------------------------------------------------------- 2. GC content and length of each group
fa <- readDNAStringSet(file.path(NF, "results_pilot/genome/genome.transcripts.fa"))
tx <- data.table(transcript_id = sub(" .*", "", names(fa)), gc = letterFrequency(fa, "GC", as.prob = TRUE)[, 1], len = width(fa))
rm(fa); invisible(gc())
gg <- merge(tx, fread(file.path(ORTHO, "tx2group.tsv")), by = "transcript_id")[, .(gc = mean(gc), len = median(len)), by = group_id]
fwrite(gg, file.path(OUT, "group_gc_length.tsv.gz"), sep = "\t")

bias_contrasts <- data.table(
  contrast = c("T_vs_W", "M_vs_W", "T_vs_M", "Mt_D1_vs_Mt_W", "Mt_D3_vs_Mt_W", "Mm_W_vs_Mt_W", "Mt_D3_vs_Mt_D1", "Fw_D1_vs_Fw_W"),
  runs = c("crosses runs (45 vs 44)", "crosses runs (45 vs 44)", "one run (45)", "crosses runs (51/52 vs 61/63/73)",
           "crosses runs (51/52 vs 61/63/73)", "crosses runs (51/52 vs 61/63/73)", "one run set (51/52)", "one run set (51/52)"),
  tissue = c("bud", "bud", "bud", "leaf", "leaf", "leaf", "leaf", "leaf"))
bias_contrasts[, crosses := grepl("^crosses", runs)]
bias <- list(); bins <- list()
for (g in c("8008", "5000")) for (i in seq_len(nrow(bias_contrasts))) {
  k <- bias_contrasts$contrast[i]
  x <- merge(rd_de(g, k)[baseMean >= 20 & !is.na(log2FC_shrunk)], gg, by.x = "id", by.y = "group_id")
  bias[[length(bias) + 1]] <- data.table(genotype = g, bias_contrasts[i], n = nrow(x),
    rho_gc = round(cor(x$log2FC_shrunk, x$gc, method = "spearman"), 3),
    rho_length = round(cor(x$log2FC_shrunk, log10(x$len), method = "spearman"), 3),
    sd_lfc = round(sd(x$log2FC_shrunk), 3))
  x[, gc_bin := cut(gc, unique(quantile(gc, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
  bins[[length(bins) + 1]] <- x[, .(gc = median(gc), median_lfc = median(log2FC_shrunk), n = .N), by = gc_bin][
    , `:=`(genotype = g, contrast = k, runs = bias_contrasts$runs[i], crosses = bias_contrasts$crosses[i])]
}
bias <- rbindlist(bias); bins <- rbindlist(bins)
fwrite(bias, file.path(OUT, "gc_length_bias.tsv"), sep = "\t"); print(bias)
bins[, contrast := factor(contrast, bias_contrasts$contrast)]
rc <- c("TRUE" = P$dir[["up"]], "FALSE" = P$dir[["down"]])
p2 <- ggplot(bins, aes(100 * gc, median_lfc, colour = as.character(crosses), group = contrast)) +
  geom_hline(yintercept = 0, colour = P$muted, linewidth = 0.3) + geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  facet_grid(genotype ~ contrast) +
  scale_colour_manual(values = rc, labels = c("TRUE" = "crosses sequencing runs", "FALSE" = "inside one run set"), name = NULL) +
  labs(x = "GC content of the group (%; 10 bins)", y = "Median shrunk log2FC",
       title = "Does the change between groups depend on GC content?",
       subtitle = "A run effect often shows as a GC trend. Groups with baseMean >= 20.") +
  theme(legend.position = "top", axis.text.x = element_text(size = 6))
save_plot(p2, file.path(OUT, "gc_bias.png"), 14, 5.5)

# ---------------------------------------------------------------- 3. genes that should be stable
ref_agi <- c(PP2AA3 = "AT1G13320", SAND = "AT2G28390", TIP41 = "AT4G34270", UBC21 = "AT5G25760", YLS8 = "AT5G08290",
             AP2M = "AT5G46630", "F-box" = "AT5G15710", UBC9 = "AT4G27960", ACT2 = "AT3G18780", ACT7 = "AT5G09810",
             UBQ10 = "AT4G05320", EF1A = "AT5G60390", GAPC1 = "AT3G04120", GAPC2 = "AT1G13440", TUB4 = "AT5G44340",
             CUL1 = "AT4G02570", EIF4A1 = "AT3G13920")
gi <- fread(file.path(ORTHO, "group_info.tsv"))
gi_l <- gi[Athaliana != "", .(agi = unlist(strsplit(Athaliana, ","))), by = group_id]
ref <- merge(gi_l, data.table(gene = names(ref_agi), agi = ref_agi), by = "agi")
ctrl <- list()
for (g in c("8008", "5000")) {
  nl <- fread(file.path(AN, sprintf("results/de/group/normcounts_%s_leaf.tsv.gz", g)))
  nb <- fread(file.path(AN, sprintf("results/de/group/normcounts_%s_bud.tsv.gz", g)))
  ml <- as.matrix(nl[, -1]); mb <- as.matrix(nb[, -1])
  st <- data.table(group_id = nl$id, leaf_mean = rowMeans(ml), leaf_cv = matrixStats::rowSds(ml) / rowMeans(ml))
  st <- merge(st, data.table(group_id = nb$id, bud_mean = rowMeans(mb)), by = "group_id")
  emp <- st[leaf_mean >= 50 & bud_mean >= 50][order(leaf_cv)][1:1000]
  ctrl[[g]] <- rbind(data.table(genotype = g, type = "empirical (least variable in leaves)", emp[, .(group_id, gene = "", leaf_cv)]),
                     data.table(genotype = g, type = "reference genes (literature)",
                                merge(ref[, .(group_id, gene)], st[, .(group_id, leaf_cv)], by = "group_id")))
}
ctrl <- rbindlist(ctrl); fwrite(ctrl, file.path(OUT, "control_genes.tsv"), sep = "\t")

ctl_lfc <- list()
for (g in c("8008", "5000")) for (k in c("T_vs_W", "M_vs_W", "T_vs_M", "Mt_D3_vs_Mt_D1", "Fw_D1_vs_Fw_W")) {
  x <- rd_de(g, k)[!is.na(padj)]
  for (ty in unique(ctrl$type)) {
    s <- x[id %in% ctrl[genotype == g & type == ty, group_id]]
    ctl_lfc[[length(ctl_lfc) + 1]] <- data.table(genotype = g, contrast = k, genes = ty, id = s$id, log2FC = s$log2FC, padj = s$padj)
  }
  ctl_lfc[[length(ctl_lfc) + 1]] <- data.table(genotype = g, contrast = k, genes = "all tested groups", id = x$id, log2FC = x$log2FC, padj = x$padj)
}
ctl_lfc <- rbindlist(ctl_lfc)
ctl_sum <- ctl_lfc[, .(n = .N, median_lfc = round(median(log2FC), 3), sd_lfc = round(sd(log2FC), 3),
                       pct_abs_lfc_gt_0.5 = round(100 * mean(abs(log2FC) > 0.5), 1),
                       pct_padj_lt_0.05 = round(100 * mean(padj < PADJ), 1)), by = .(genotype, contrast, genes)]
fwrite(ctl_sum, file.path(OUT, "control_genes_lfc_summary.tsv"), sep = "\t"); print(ctl_sum)
refrows <- ctl_lfc[genes == "reference genes (literature)"]
refrows <- merge(refrows, unique(ref[, .(id = group_id, gene)]), by = "id", allow.cartesian = TRUE)
fwrite(refrows, file.path(OUT, "reference_genes_lfc.tsv"), sep = "\t")
cl <- ctl_lfc[genes != "all tested groups" & contrast %in% c("T_vs_W", "M_vs_W", "T_vs_M")]
cl[, contrast := factor(contrast, c("T_vs_W", "M_vs_W", "T_vs_M"),
                        c("T vs W (crosses runs)", "M vs W (crosses runs)", "T vs M (one run)"))]
p3 <- ggplot(cl, aes(contrast, pmax(pmin(log2FC, 3), -3), fill = contrast)) +
  geom_hline(yintercept = 0, colour = P$muted, linewidth = 0.3) +
  geom_violin(colour = NA, alpha = 0.6, scale = "width") + geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_grid(genes ~ genotype, scales = "free_y") +
  scale_fill_manual(values = setNames(c(P$dir[["up"]], P$dir[["up"]], P$dir[["down"]]), levels(cl$contrast)), guide = "none") +
  labs(x = NULL, y = "log2FC (clamped to ±3)", title = "Do genes that should be stable move between bud groups?",
       subtitle = "If the run effect were small, the two run-crossing contrasts would look like T vs M.")
save_plot(p3, file.path(OUT, "control_genes_lfc.png"), 9, 6.5)

# ---------------------------------------------------------------- 4. RUVg and 5. leaf drought programme in buds
txi <- readRDS(file.path(AN, "data/txi_group.rds"))
ruv_sum <- list(); ruv_w <- list(); prog <- list()
for (g in c("8008", "5000")) {
  md <- meta[meta$set == "bud" & meta$genotype == g, ]; md$condition <- factor(md$condition, c("W", "T", "M"))
  tx_s <- subset_txi(txi, md$sample)
  ctl_g <- intersect(ctrl[genotype == g & grepl("^empirical", type), group_id], rownames(tx_s$counts))
  es <- newSeqExpressionSet(round(tx_s$counts), phenoData = AnnotatedDataFrame(data.frame(condition = md$condition, row.names = md$sample)))
  es <- betweenLaneNormalization(es, which = "upper")
  leaf_d1 <- rd_de(g, "Mt_D1_vs_Mt_W")
  lsets <- list(leaf_induced = leaf_d1[deg == "up", id], leaf_repressed = leaf_d1[deg == "down", id])
  for (k in 0:2) {
    cd <- md
    if (k > 0) { r <- RUVg(es, ctl_g, k = k); for (j in seq_len(k)) cd[[paste0("W_", j)]] <- pData(r)[[paste0("W_", j)]] }
    form <- as.formula(paste("~", paste(c(if (k > 0) paste0("W_", seq_len(k)), "condition"), collapse = " + ")))
    dd <- DESeq(DESeqDataSetFromTximport(tx_s, cd, form), quiet = TRUE)
    for (cc in list(c("T", "W"), c("M", "W"), c("T", "M"))) {
      rr <- results(dd, contrast = c("condition", cc), alpha = PADJ)
      dg <- call_deg(rr$padj, rr$log2FoldChange); nm <- paste(cc, collapse = "_vs_")
      ruv_sum[[length(ruv_sum) + 1]] <- data.table(genotype = g, k = k, contrast = nm, up = sum(dg == "up"), down = sum(dg == "down"))
      st <- setNames(rr$stat, rownames(rr)); st <- st[!is.na(st)]
      cp <- cameraPR(st, lapply(lsets, function(s) which(names(st) %in% s)))
      prog[[length(prog) + 1]] <- data.table(genotype = g, k = k, contrast = nm, set = rownames(cp), cp)
    }
    if (k > 0) ruv_w[[length(ruv_w) + 1]] <- data.table(genotype = g, k = k, sample = md$sample, condition = md$condition,
                                                         run = md$run, cd[, paste0("W_", seq_len(k)), drop = FALSE])
  }
}
ruv_sum <- rbindlist(ruv_sum); fwrite(ruv_sum, file.path(OUT, "ruvg_deg_summary.tsv"), sep = "\t"); print(ruv_sum)
ruv_w <- rbindlist(ruv_w, fill = TRUE); fwrite(ruv_w, file.path(OUT, "ruvg_factors.tsv"), sep = "\t")
prog <- rbindlist(prog); fwrite(prog, file.path(OUT, "leaf_drought_programme_in_buds.tsv"), sep = "\t"); print(prog)
wc <- ruv_w[k == 1, .(cor_W1_with_watered = round(cor(W_1, as.numeric(condition == "W")), 3)), by = genotype]
fwrite(wc, file.path(OUT, "ruvg_W1_vs_run.tsv"), sep = "\t"); print(wc)
# correlations between bud and leaf drought contrasts (Wald statistic, Spearman)
cors <- rbindlist(lapply(c("8008", "5000"), function(g) {
  ks <- c("T_vs_W", "M_vs_W", "T_vs_M", "Mt_D1_vs_Mt_W")
  m <- Reduce(function(a, b) merge(a, b, by = "id"), lapply(ks, function(k) setnames(rd_de(g, k)[!is.na(stat), .(id, stat)], "stat", k)))
  cc <- cor(as.matrix(m[, -1]), method = "spearman")
  data.table(genotype = g, a = rep(ks, each = 4), b = rep(ks, 4), rho = round(as.vector(cc), 3))[a < b]
}))
fwrite(cors, file.path(OUT, "bud_leaf_correlations.tsv"), sep = "\t"); print(cors)
rs <- copy(ruv_sum)[, contrast := factor(contrast, c("T_vs_W", "M_vs_W", "T_vs_M"))]
rs <- melt(rs, id.vars = c("genotype", "k", "contrast"), variable.name = "direction", value.name = "DEGs")
p4 <- ggplot(rs, aes(factor(k), DEGs, fill = direction)) + geom_col(width = 0.7) +
  facet_grid(genotype ~ contrast) + scale_fill_manual(values = P$dir[c("up", "down")], name = NULL) +
  labs(x = "Unwanted factors removed by RUVg (k)", y = "DEGs", title = "Bud DEGs before (k = 0) and after RUVg",
       subtitle = "Controls: 1,000 groups that vary least across the genotype's leaf libraries")
save_plot(p4, file.path(OUT, "ruvg_degs.png"), 8, 5.5)
cat("done\n")
