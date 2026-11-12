#!/usr/bin/env Rscript
# Genotype identity check from expression. POJ2878 is allele-resolved, and the two cultivars carry different
# alleles, so genes expressed in some libraries and silent in others (presence/absence) mark genotype regardless
# of tissue or treatment. Libraries are clustered on these on/off patterns and compared with their labels.
suppressMessages({ library(data.table); library(pheatmap) })
D <- "/Storage/data1/jorge.munoz/DOLORES/analysis"; OUT <- file.path(D, "results/qc/genotype_check")
meta <- fread(file.path(D, "data/samples.tsv"), colClasses = list(character = "genotype")); setDF(meta); rownames(meta) <- meta$sample
txi <- readRDS(file.path(D, "data/txi_gene.rds")); cnt <- txi$counts[, meta$sample]
cpm <- t(t(cnt) / colSums(cnt)) * 1e6
on <- cpm >= 5; off <- cnt == 0
# informative: clearly on in >= 8 libraries and fully off in >= 8 libraries
inf <- rowSums(on) >= 8 & rowSums(off) >= 8
cat("presence/absence genes:", sum(inf), "\n")
b <- (cpm[inf, ] > 1) * 1
d <- dist(t(b), method = "binary"); hc <- hclust(d, "average")
k2 <- cutree(hc, 2)
tab <- table(cluster = k2, label = meta$genotype); print(tab)
maj <- apply(tab, 1, function(r) colnames(tab)[which.max(r)])
meta$pa_cluster_genotype <- maj[as.character(k2)]
# score: correlation of each library's on/off profile with each genotype's consensus (leave-one-out)
score <- t(sapply(seq_len(ncol(b)), function(i) {
  sapply(c("5000", "8008"), function(g) { j <- setdiff(which(meta$genotype == g), i); cor(b[, i], rowMeans(b[, j])) })
}))
meta$cor_5000 <- round(score[, 1], 3); meta$cor_8008 <- round(score[, 2], 3)
meta$best_match <- ifelse(meta$cor_5000 > meta$cor_8008, "5000", "8008")
meta$mismatch <- meta$best_match != meta$genotype
fwrite(meta[, c("sample", "set", "condition", "genotype", "pa_cluster_genotype", "cor_5000", "cor_8008", "best_match", "mismatch")],
       file.path(OUT, "genotype_check.tsv"), sep = "\t")
print(meta[meta$mismatch | meta$pa_cluster_genotype != meta$genotype, c("sample", "genotype", "pa_cluster_genotype", "cor_5000", "cor_8008")])
pheatmap(as.matrix(1 - d), annotation_col = meta[, c("set", "genotype")], fontsize = 6,
         clustering_distance_rows = d, clustering_distance_cols = d, filename = file.path(OUT, "presence_absence_similarity.png"), width = 13, height = 11)
