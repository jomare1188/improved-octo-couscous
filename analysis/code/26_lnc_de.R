#!/usr/bin/env Rscript
# Genome-wide lncRNA search, step 5 of 5: do lncRNAs respond to drought, and do any behave as memory genes?
# Salmon counts on annotated transcripts + lncRNA candidates (25_lnc_quant.sh). Coding transcripts are summed
# into their allele groups and lncRNA candidates into their clusters of allele copies (24_lnc_candidates.sh),
# so the lncRNAs are tested among all genes, with the same models, contrasts and thresholds as 03_de.R.
# Output: results/lncrna_genome/ (only the lncRNA rows).   qsub -N lncde -hold_jid lncquant code/run_r.sh code/26_lnc_de.R
source("code/common.R")
suppressMessages({ library(DESeq2); library(tximport); library(BiocParallel) }); register(SerialParam())
W   <- "/Storage/data2/jorge.munoz/lnc_genome"
OUT <- file.path(AN, "results/lncrna_genome"); dir.create(OUT, showWarnings = FALSE)
meta <- read_meta()
files <- file.path(W, "quant", meta$library, "quant.sf"); names(files) <- meta$sample
stopifnot(all(file.exists(files)))
t2g <- rbind(fread(file.path(ORTHO, "tx2group.tsv"))[, .(transcript_id, group_id)], fread(file.path(W, "lnc_tx2group.tsv")))
txi <- tximport(files, type = "salmon", tx2gene = as.data.frame(t2g), countsFromAbundance = "no")
saveRDS(txi, file.path(W, "txi_group_lnc.rds"))
is_lnc <- function(x) grepl("^LNC:", x)

# one row per lncRNA cluster
cand <- fread(file.path(W, "candidates.tsv"))[lncRNA == 1]
cl <- merge(fread(file.path(W, "lnc_tx2group.tsv")), cand, by.x = "transcript_id", by.y = "transcript")
info <- cl[, .(n_copies = .N, class_code = names(sort(table(class_code), decreasing = TRUE))[1], max_exons = max(exons),
               length_representative = length[match(sub("^LNC:", "", group_id[1]), transcript_id)]), by = group_id]
ab <- txi$abundance[is_lnc(rownames(txi$abundance)), ]
info <- merge(info, data.table(group_id = rownames(ab), libs_tpm_ge1 = rowSums(ab >= 1), max_tpm = round(apply(ab, 1, max), 2)), by = "group_id", all.x = TRUE)
info[, high_confidence := max_exons >= 2 & libs_tpm_ge1 >= 2]

contrasts_leaf <- list(c("Mt_D1", "Mt_W"), c("Mt_D3", "Mt_W"), c("Mt_D3", "Mt_D1"),
                       c("Fw_D1", "Fw_W"), c("Ft_D1", "Fw_D1"), c("Fm_D1", "Fw_D1"))
contrasts_bud <- list(c("T", "M"), c("T", "W"), c("M", "W"))
summ <- list(); allres <- list(); mem <- list()
for (g in c("8008", "5000")) for (st in c("leaf", "bud")) {
  m <- meta[meta$genotype == g & (if (st == "bud") meta$set == "bud" else meta$set != "bud"), ]
  m$condition <- factor(m$condition)
  dds <- DESeq(DESeqDataSetFromTximport(subset_txi(txi, m$sample), m, ~ condition), quiet = TRUE)
  ns <- table(m$condition); rl <- list()
  for (cc in if (st == "bud") contrasts_bud else contrasts_leaf) {
    if (!all(cc %in% names(ns))) next
    r <- results(dds, contrast = c("condition", cc), alpha = PADJ)     # padj over all genes, coding and lncRNA
    name <- sprintf("%s_vs_%s", cc[1], cc[2])
    dt <- data.table(genotype = g, contrast = name, id = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, padj = r$padj)
    dt[, deg := call_deg(padj, log2FC)]
    status <- if (st == "bud" && "W" %in% cc) "confounded with run" else if (min(ns[cc]) < 2) "exploratory" else "conclusive"
    summ[[length(summ) + 1]] <- data.table(genotype = g, contrast = name, n1 = ns[[cc[1]]], n2 = ns[[cc[2]]], status = status,
      lnc_tested = sum(is_lnc(dt$id) & !is.na(dt$padj)), lnc_up = sum(is_lnc(dt$id) & dt$deg == "up"), lnc_down = sum(is_lnc(dt$id) & dt$deg == "down"),
      lnc_hc_up = sum(dt$id %in% info[high_confidence == TRUE, group_id] & dt$deg == "up"),
      lnc_hc_down = sum(dt$id %in% info[high_confidence == TRUE, group_id] & dt$deg == "down"),
      coding_up = sum(!is_lnc(dt$id) & dt$deg == "up"), coding_down = sum(!is_lnc(dt$id) & dt$deg == "down"))
    rl[[name]] <- dt; allres[[length(allres) + 1]] <- dt[is_lnc(id)]
  }
  if (st == "leaf" && all(c("Mt_D1_vs_Mt_W", "Mt_D3_vs_Mt_D1", "Mt_D3_vs_Mt_W") %in% names(rl))) {
    a <- rl[["Mt_D1_vs_Mt_W"]]; b <- rl[["Mt_D3_vs_Mt_D1"]]; d3 <- rl[["Mt_D3_vs_Mt_W"]]
    sgn <- function(x) fifelse(x$deg == "up", "+", fifelse(x$deg == "down", "-", "="))
    x <- data.table(genotype = g, id = a$id, D1_vs_W = a$log2FC, D1_vs_W_padj = a$padj, D3_vs_D1 = b$log2FC, D3_vs_D1_padj = b$padj,
                    D3_vs_W = d3$log2FC, D3_vs_W_padj = d3$padj, r1 = sgn(a), r2 = sgn(b), r3 = sgn(d3))
    x[, class := fifelse(r1 != "=", paste0(r1, r2), fifelse(r3 != "=", paste0("0", r3), "none"))]
    mem[[g]] <- x[is_lnc(id), !c("r1", "r2", "r3")]
  }
}
summ <- rbindlist(summ); allres <- rbindlist(allres); mem <- rbindlist(mem)
mem <- merge(mem, info, by.x = "id", by.y = "group_id", all.x = TRUE)
fwrite(info, file.path(OUT, "lncRNA_clusters.tsv"), sep = "\t")
fwrite(summ, file.path(OUT, "lnc_deg_summary.tsv"), sep = "\t")
fwrite(allres, file.path(OUT, "lnc_de_all_contrasts.tsv.gz"), sep = "\t")
fwrite(mem, file.path(OUT, "lnc_memory_classes_tillering.tsv.gz"), sep = "\t")
cat("lncRNA clusters:", nrow(info), " high confidence:", sum(info$high_confidence, na.rm = TRUE), "\n")
print(summ)
print(dcast(mem[, .N, by = .(genotype, class, hc = high_confidence)], class ~ genotype + hc, value.var = "N", fill = 0))
print(mem[genotype == "8008" & class %in% c("++", "--", "+-", "-+")][order(class, -abs(D3_vs_D1)), .(id, class, D1_vs_W = round(D1_vs_W, 2), D3_vs_D1 = round(D3_vs_D1, 2), n_copies, class_code, max_exons, max_tpm)])
