#!/usr/bin/env Rscript
# Input for WGCNA (co-expression modules), allele-group level.
# Leaves only (32 libraries, both genotypes): buds are too few (16) and their watered group is run-confounded.
#  - VST of counts (groups with >= 10 reads in >= 3 leaf libraries and mean normalised count >= 10)
#  - the genotype main effect is removed (limma removeBatchEffect, keeping the condition means), so modules
#    capture co-variation shared by both genotypes rather than "5000 vs 8008" differences
#  - the 8,000 most variable groups after that correction
#  - traits per library: current drought, repeated drought (D3), maternal stress (Ft/Fm), propagule, maturation,
#    late run (Mt_W, the only group on runs 61/63/73)
#  - bud VST of the same groups (centred per genotype), to project the leaf modules onto buds later
# Run from analysis/:  qsub -N wgprep code/run_r.sh code/14_wgcna_prep.R
source("code/common.R")
suppressMessages({ library(DESeq2); library(limma); library(matrixStats) })
OUT <- file.path(AN, "results/wgcna"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NTOP <- 8000
meta <- read_meta(); txi <- readRDS(file.path(AN, "data/txi_group.rds"))

leaf <- meta[meta$set != "bud", ]
dds <- DESeqDataSetFromTximport(subset_txi(txi, leaf$sample, min_n = 3), leaf, ~ 1)
dds <- estimateSizeFactors(dds)
dds <- dds[rowMeans(counts(dds, normalized = TRUE)) >= 10, ]
x <- assay(vst(dds, blind = TRUE))
cond <- factor(leaf$condition)
xc <- removeBatchEffect(x, batch = leaf$genotype, design = model.matrix(~ cond))
top <- order(rowVars(xc), decreasing = TRUE)[seq_len(NTOP)]
expr <- t(xc[top, ])
saveRDS(expr, file.path(OUT, "input_expr.rds"))
fwrite(data.table(group_id = colnames(expr), variance = rowVars(xc)[top]), file.path(OUT, "input_groups.tsv"), sep = "\t")
cat("leaf groups after filter:", nrow(x), " used:", ncol(expr), " libraries:", nrow(expr), "\n")

tr <- data.table(sample = leaf$sample, genotype = leaf$genotype, condition = leaf$condition, run = leaf$run)
tr[, `:=`(current   = as.integer(grepl("_D[13]$", condition)),
          repeated  = as.integer(grepl("_D3$", condition)),
          maternal  = as.integer(grepl("^F[tm]_", condition)),
          propagule = as.integer(grepl("^F", condition)),
          maturation = as.integer(grepl("^Mm_", condition)),
          late_run  = as.integer(condition == "Mt_W"),
          g8008     = as.integer(genotype == "8008"))]
fwrite(tr, file.path(OUT, "traits.tsv"), sep = "\t"); print(tr[, .N, by = .(condition, current, repeated, maternal, propagule, maturation, late_run)])

bud <- meta[meta$set == "bud", ]
db <- DESeqDataSetFromTximport(subset_txi(txi, bud$sample, min_n = 2), bud, ~ 1)
xb <- assay(vst(db, blind = TRUE))
xb <- xb[intersect(colnames(expr), rownames(xb)), ]
for (g in unique(bud$genotype)) { i <- bud$genotype == g; xb[, i] <- xb[, i] - rowMeans(xb[, i]) }
saveRDS(xb, file.path(OUT, "bud_expr_centred.rds"))
cat("bud groups for projection:", nrow(xb), "\n")
