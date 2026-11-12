#!/usr/bin/env Rscript
# Differential expression with DESeq2, each genotype analysed separately.   Rscript 03_de.R <gene|group>
#  - buds:   ~ condition (W, T = tillering-stressed, M = maturation-stressed), n = 3 per group
#  - leaves: mother and propagule leaves in one model, ~ condition, so groups with one replicate borrow
#            dispersion from the others. Contrasts are pulled per hypothesis.
# Every contrast table carries the group sizes and a status:
#   conclusive  = both groups have >= 2 replicates
#   exploratory = a group has 1 replicate (informs, no conclusions)
#   descriptive = maturation mother leaves (agreed: descriptive only)
# DEG = padj < 0.05 and |log2FC| >= 1 (unshrunk LFC; a shrunk LFC column, type "normal", is given for ranking).
# Memory classes (tillering mother leaves; Ding et al. 2013/2014 definition):
#   first sign = response at D1 vs W, second sign = change D3 vs D1; '=' = no significant change.
#   Late responders (no D1 response, D3 vs W significant) are classed 0+ / 0-.
suppressMessages({ library(DESeq2); library(data.table); library(BiocParallel) })
register(SerialParam())   # forked workers hung at group level; serial is fast enough
level <- commandArgs(TRUE)[1]; if (is.na(level)) level <- "gene"
D   <- "/Storage/data1/jorge.munoz/DOLORES/analysis"
OUT <- file.path(D, "results/de", level); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta <- fread(file.path(D, "data/samples.tsv"), colClasses = list(character = "genotype")); setDF(meta); rownames(meta) <- meta$sample
txi  <- readRDS(file.path(D, sprintf("data/txi_%s.rds", level)))
# libraries whose genotype (from expressed alleles) does not match the label: left out until resolved
excl <- fread(file.path(D, "data/excluded_samples.tsv"))$sample
meta <- meta[!meta$sample %in% excl, ]
cat("excluded:", excl, "\n")
PADJ <- 0.05; LFC <- 1

contrasts_leaf <- list(
  # H1: memory across cycles, tillering mother leaves
  c("Mt_D1", "Mt_W"), c("Mt_D3", "Mt_W"), c("Mt_D3", "Mt_D1"),
  # maturation mother leaves (descriptive)
  c("Mm_D1", "Mm_W"), c("Mm_D3", "Mm_W"), c("Mm_D3", "Mm_D1"), c("Mm_W", "Mt_W"),
  # H4: propagules
  c("Fw_D1", "Fw_W"), c("Ft_D1", "Fw_D1"), c("Fm_D1", "Fw_D1"), c("Ft_D1", "Fw_W"), c("Fm_D1", "Fw_W"))
contrasts_bud <- list(c("T", "W"), c("M", "W"), c("T", "M"))

subset_txi <- function(idx) {   # idx = sample names
  k <- rowSums(txi$counts[, idx] >= 10) >= 2
  list(abundance = txi$abundance[k, idx], counts = txi$counts[k, idx], length = txi$length[k, idx],
       countsFromAbundance = "no")
}

summary_rows <- list()
run_set <- function(set_name, sets, genotype, cons) {
  idx <- which(meta$set %in% sets & meta$genotype == genotype)
  md  <- meta[idx, ]; md$condition <- factor(md$condition)
  dds <- DESeqDataSetFromTximport(subset_txi(md$sample), md, ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  tag <- sprintf("%s_%s", genotype, set_name)
  saveRDS(dds, file.path(OUT, sprintf("dds_%s.rds", tag)))
  fwrite(data.table(id = rownames(dds), counts(dds, normalized = TRUE)), file.path(OUT, sprintf("normcounts_%s.tsv.gz", tag)), sep = "\t")
  ns <- table(md$condition)
  res_list <- list()
  for (cc in cons) {
    if (!all(cc %in% names(ns))) next
    n1 <- ns[[cc[1]]]; n2 <- ns[[cc[2]]]
    status <- if (grepl("^Mm", cc[1]) || grepl("^Mm", cc[2])) "descriptive" else if (min(n1, n2) < 2) "exploratory" else "conclusive"
    r  <- results(dds, contrast = c("condition", cc), alpha = PADJ)
    sh <- lfcShrink(dds, contrast = c("condition", cc), res = r, type = "normal", quiet = TRUE)
    name <- sprintf("%s_vs_%s", cc[1], cc[2])
    dt <- data.table(id = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                     log2FC_shrunk = sh$log2FoldChange, stat = r$stat, pvalue = r$pvalue, padj = r$padj)
    dt[, deg := fifelse(!is.na(padj) & padj < PADJ & abs(log2FC) >= LFC, fifelse(log2FC > 0, "up", "down"), "")]
    attr(dt, "status") <- status
    fwrite(dt, file.path(OUT, sprintf("de_%s_%s.tsv.gz", genotype, name)), sep = "\t")
    res_list[[name]] <- dt
    summary_rows[[length(summary_rows) + 1]] <<- data.table(genotype = genotype, set = set_name, contrast = name,
      n1 = n1, n2 = n2, status = status, tested = sum(!is.na(dt$padj)), up = sum(dt$deg == "up"), down = sum(dt$deg == "down"))
  }
  res_list
}

for (g in c("5000", "8008")) {
  run_set("bud", "bud", g, contrasts_bud)
  rl <- run_set("leaf", c("mother_leaf", "propagule_leaf"), g, contrasts_leaf)
  # memory classes, tillering mother leaves
  a <- rl[["Mt_D1_vs_Mt_W"]]; b <- rl[["Mt_D3_vs_Mt_D1"]]; d3 <- rl[["Mt_D3_vs_Mt_W"]]
  if (!is.null(a) && !is.null(b)) {
    sgn <- function(dt) fifelse(dt$deg == "up", "+", fifelse(dt$deg == "down", "-", "="))
    m <- data.table(id = a$id, D1_vs_W = a$log2FC, D1_vs_W_padj = a$padj, D3_vs_D1 = b$log2FC, D3_vs_D1_padj = b$padj,
                    D3_vs_W = d3$log2FC, D3_vs_W_padj = d3$padj, r1 = sgn(a), r2 = sgn(b), r3 = sgn(d3))
    m[, class := fifelse(r1 != "=", paste0(r1, r2), fifelse(r3 != "=", paste0("0", r3), "none"))]
    m[, status := fifelse(g == "5000", "exploratory (Mt_D1 n=1)", "conclusive")]
    fwrite(m[, !c("r1", "r2", "r3")], file.path(OUT, sprintf("memory_classes_%s_tillering.tsv.gz", g)), sep = "\t")
    s <- m[, .N, by = class][order(class)]; s[, genotype := g]
    fwrite(s, file.path(OUT, sprintf("memory_class_counts_%s_tillering.tsv", g)), sep = "\t")
    print(s)
  }
}
sm <- rbindlist(summary_rows)
fwrite(sm, file.path(OUT, "deg_summary.tsv"), sep = "\t")
print(sm)
