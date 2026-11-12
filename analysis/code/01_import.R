#!/usr/bin/env Rscript
# Import Salmon quantifications (51 libraries) with tximport and sum transcripts into features.
#   Rscript 01_import.R <level>      level = gene   (POJ2878 gene = one allele; Salmon tx2gene)
#                                    level = group  (allele groups anchored on sorghum, from orthology/)
# Output: data/txi_<level>.rds (tximport list: counts, abundance, length) and data/samples.tsv
suppressMessages({ library(tximport); library(data.table) })
level <- commandArgs(TRUE)[1]; if (is.na(level)) level <- "gene"
D   <- "/Storage/data1/jorge.munoz/DOLORES"
NF  <- file.path(D, "nf-rnaseq/run_POJ2878")
OUT <- file.path(D, "analysis/data"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

meta <- fread(file.path(NF, "sample_metadata_POJ2878_clean.csv"))
meta[, condition := sub("^[0-9]+_", "", group)]                 # e.g. Mt_D1, Tillering_D
meta[set == "Mother bud", condition := c(Maturation_W = "W", Tillering_D = "T", Maturation_D = "M")[condition]]
meta[, set := c("Mother bud" = "bud", "Mother leaf" = "mother_leaf", "Propagule leaf" = "propagule_leaf")[set]]
meta[, genotype := as.character(genotype)]
meta[, n_in_group := .N, by = group]
fwrite(meta, file.path(OUT, "samples.tsv"), sep = "\t")

files <- file.path(NF, "results/salmon", meta$sample, "quant.sf"); names(files) <- meta$sample
stopifnot(all(file.exists(files)))

if (level == "gene") {
  t2g <- fread(file.path(NF, "results/salmon/salmon.merged.tx2gene.tsv"))[, .(transcript_id, gene_id)]
} else if (level == "group") {
  t2g <- fread(file.path(D, "orthology/allele_groups/tx2group.tsv"))[, .(transcript_id, group_id)]
} else stop("level must be gene or group")

txi <- tximport(files, type = "salmon", tx2gene = as.data.frame(t2g), countsFromAbundance = "no")
saveRDS(txi, file.path(OUT, sprintf("txi_%s.rds", level)))
cat(level, ": ", nrow(txi$counts), " features x ", ncol(txi$counts), " samples\n", sep = "")
cat("features with >= 10 reads in >= 2 samples: ", sum(rowSums(txi$counts >= 10) >= 2), "\n")
