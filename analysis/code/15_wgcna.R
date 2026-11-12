#!/usr/bin/env Rscript
# WGCNA network and modules on the leaf input of 14_wgcna_prep.R (signed network, biweight midcorrelation).
# Needs the r_wgcna env (WGCNA is not in r_RNASEQ):
#   qsub -N wgcna -v RENV=/Storage/data2/jorge.munoz/conda_envs/r_wgcna code/run_r.sh code/15_wgcna.R
# Soft-threshold power: the lowest power with signed scale-free fit R^2 >= 0.80; if none reaches it, 14
# (WGCNA FAQ value for signed networks with 30-40 samples).
suppressMessages({ library(WGCNA); library(data.table) })
options(stringsAsFactors = FALSE)
enableWGCNAThreads(4)
OUT <- "/Storage/data1/jorge.munoz/DOLORES/analysis/results/wgcna"
expr <- readRDS(file.path(OUT, "input_expr.rds"))
gsg <- goodSamplesGenes(expr, verbose = 0)
if (!gsg$allOK) expr <- expr[gsg$goodSamples, gsg$goodGenes]
cat("libraries:", nrow(expr), " groups:", ncol(expr), "\n")

png(file.path(OUT, "sample_tree.png"), width = 2400, height = 1200, res = 200)
plot(hclust(dist(expr), method = "average"), main = "Leaf libraries clustered on the WGCNA input", sub = "", xlab = "", cex = 0.6)
dev.off()

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(expr, powerVector = powers, networkType = "signed", corFnc = bicor,
                         corOptions = list(maxPOutliers = 0.05), verbose = 0)
fit <- as.data.table(sft$fitIndices)[, signed_R2 := -sign(slope) * SFT.R.sq]
pw <- fit$Power[which(fit$signed_R2 >= 0.80)[1]]
rule <- "lowest power with signed R2 >= 0.80"
if (is.na(pw)) { pw <- 14; rule <- "no power reached R2 0.80; FAQ value for signed networks, 30-40 samples" }
fwrite(fit, file.path(OUT, "soft_threshold.tsv"), sep = "\t"); print(fit)
writeLines(c(paste("power", pw), paste("rule", rule)), file.path(OUT, "power.txt")); cat("power", pw, "(", rule, ")\n")

net <- blockwiseModules(expr, power = pw, networkType = "signed", TOMType = "signed", corType = "bicor",
                        maxPOutliers = 0.05, minModuleSize = 30, deepSplit = 2, mergeCutHeight = 0.25,
                        pamRespectsDendro = FALSE, numericLabels = TRUE, maxBlockSize = 10000,
                        saveTOMs = FALSE, nThreads = 4, verbose = 2)
lab <- net$colors
# rename modules M01, M02 ... by size (M00 = unassigned)
sizes <- sort(table(lab[lab != 0]), decreasing = TRUE)
newlab <- setNames(sprintf("M%02d", seq_along(sizes)), names(sizes)); newlab["0"] <- "M00"
mod <- newlab[as.character(lab)]
MEs <- moduleEigengenes(expr, colors = mod, excludeGrey = TRUE)$eigengenes
colnames(MEs) <- sub("^ME", "", colnames(MEs))
kme <- as.matrix(signedKME(expr, MEs, corFnc = "bicor", corOptions = "maxPOutliers = 0.05"))
colnames(kme) <- sub("^kME", "", colnames(kme))
own <- rep(NA_real_, ncol(expr)); a <- mod != "M00"   # kME of each group in its own module
own[a] <- kme[cbind(which(a), match(mod[a], colnames(kme)))]
fwrite(data.table(group_id = colnames(expr), module = mod, kME_own = round(own, 3)), file.path(OUT, "modules.tsv"), sep = "\t")
fwrite(data.table(sample = rownames(expr), MEs), file.path(OUT, "eigengenes.tsv"), sep = "\t")
fwrite(data.table(group_id = colnames(expr), round(kme, 3)), file.path(OUT, "kme_all.tsv.gz"), sep = "\t")
print(table(mod))

png(file.path(OUT, "module_dendrogram.png"), width = 2400, height = 1300, res = 200)
plotDendroAndColors(net$dendrograms[[1]], labels2colors(lab[net$blockGenes[[1]]]), "Module", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = sprintf("Leaf co-expression modules (power %d)", pw))
dev.off()
cat("done\n")
