#!/usr/bin/env Rscript
# Do libraries of the same group differ more across sequencing runs than within one run?
# For every group with >= 2 libraries: distance between each pair of libraries (root mean square difference of
# VST allele-group expression, groups with mean VST > 6), with the runs of the two libraries.
# Used for runs 51 vs 52 (leaves), which some contrasts cross (e.g. 5000 Mt_D1 on 51, Mt_D3 on 52; the propagule
# reference Fw_D1 on the other run from Ft_D1). Runs in seconds:  Rscript code/20_run_check.R   (from analysis/)
source("code/common.R")
suppressMessages(library(SummarizedExperiment))
OUT <- file.path(AN, "results/qc"); meta <- read_meta()
vst <- assay(readRDS(file.path(AN, "data/vst_group.rds")))[, meta$sample]
pairs <- rbindlist(lapply(split(meta, meta$group), function(m) {
  if (nrow(m) < 2) return(NULL)
  x <- vst[, m$sample]; x <- x[rowMeans(x) > 6, ]
  ij <- t(combn(nrow(m), 2))
  data.table(group = m$group[1], set = m$set[1], library_1 = m$library[ij[, 1]], library_2 = m$library[ij[, 2]],
             run_1 = m$run[ij[, 1]], run_2 = m$run[ij[, 2]],
             rms_diff = round(sqrt(colMeans((x[, ij[, 1], drop = FALSE] - x[, ij[, 2], drop = FALSE])^2)), 3))
}))
pairs[, same_run := run_1 == run_2]
fwrite(pairs, file.path(OUT, "run_pair_distances.tsv"), sep = "\t")
# summary: leaves on runs 51/52, and the watered tillering leaves on the later runs
leaf5152 <- pairs[set != "bud" & run_1 %in% c("51", "52") & run_2 %in% c("51", "52")]
sm <- rbind(leaf5152[, .(what = "leaf, runs 51/52", n_pairs = .N, median_rms = median(rms_diff),
                         min = min(rms_diff), max = max(rms_diff)), by = same_run],
            pairs[grepl("Mt_W", group), .(what = "Mt_W, runs 61/63/73", n_pairs = .N, median_rms = median(rms_diff),
                                          min = min(rms_diff), max = max(rms_diff)), by = same_run])
fwrite(sm, file.path(OUT, "run_pair_distances_summary.tsv"), sep = "\t")
print(sm); print(pairs[grepl("Mt_W", group)])
