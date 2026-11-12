#!/usr/bin/env Rscript
# Known drought-memory genes in our data (allele-group level).
#  1. Map literature genes to POJ2878 allele groups through OrthoFinder orthologues:
#       Arabidopsis memory classes (Ding 2013), maize memory classes (Ding 2014; v3 IDs -> v5 with MaizeGDB),
#       curated candidates (memory_genes/curated_candidates.tsv: ABA, LEA, osmolytes, chromatin, DNA methylation).
#  2. Gene-set tests (limma cameraPR on the DESeq2 Wald statistic) of each literature set in every contrast.
#  3. Overlap of our tillering memory classes with the literature classes (Fisher tests).
#  4. Heatmap of log2FC of curated candidate groups across contrasts.
suppressMessages({ library(data.table); library(limma); library(pheatmap) })
D  <- "/Storage/data1/jorge.munoz/DOLORES"
MG <- file.path(D, "memory_genes"); DE <- file.path(D, "analysis/results/de/group")
OUT <- file.path(D, "analysis/results/memory_genes"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

gi <- fread(file.path(D, "orthology/allele_groups/group_info.tsv"))
long <- function(col, sp) { x <- gi[get(col) != "", .(orth = unlist(strsplit(get(col), ","))), by = group_id]; x[, species := sp]; x }
o2g <- rbind(long("Athaliana", "Athaliana"), long("Zmays", "Zmays"))

ding <- fread(file.path(MG, "ding_memory_classes.tsv"))
v3v5 <- fread(file.path(MG, "literature/maizegdb/B73v3_to_B73v5.tsv"), header = FALSE, col.names = c("v3", "v5"))
v3v5 <- v3v5[, .(v5 = unlist(strsplit(v5, ","))), by = v3]
ath <- unique(ding[species == "Athaliana", .(lit_species = "Athaliana", lit_gene = gene_id, orth = gene_id, symbol, class, source)])
zm  <- merge(unique(ding[species == "Zmays_v3", .(lit_gene = gene_id, class, source)]), v3v5, by.x = "lit_gene", by.y = "v3")
zm  <- zm[, .(lit_species = "Zmays", lit_gene, orth = v5, symbol = NA_character_, class, source)]
cur <- fread(file.path(MG, "curated_candidates.tsv"))
cur <- cur[, .(lit_species = "Athaliana", lit_gene = gene_id, orth = gene_id, symbol, class = paste0("curated:", category), source = why)]
lit <- rbind(ath, zm, cur)
map <- merge(lit, o2g, by.x = c("orth", "lit_species"), by.y = c("orth", "species"), allow.cartesian = TRUE)
fwrite(map, file.path(OUT, "literature_to_groups.tsv"), sep = "\t")
cat("literature genes mapped to groups:\n"); print(map[, .(genes = uniqueN(lit_gene), groups = uniqueN(group_id)), by = .(lit_species, class)][order(lit_species, class)])

# gene sets: species x memory class, and curated categories
map[, set := fifelse(grepl("^curated", class), class, paste(lit_species, class))]
sets <- split(map$group_id, map$set); sets <- lapply(sets, unique)

files <- list.files(DE, "^de_.*\\.tsv\\.gz$", full.names = TRUE)
sm <- fread(file.path(DE, "deg_summary.tsv"))
cam <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)[!is.na(stat)]
  nm <- sub("^de_", "", sub("\\.tsv\\.gz$", "", basename(f)))
  idx <- lapply(sets, function(s) which(dt$id %in% s)); idx <- idx[lengths(idx) >= 3]
  r <- cameraPR(setNames(dt$stat, dt$id), idx)
  data.table(contrast = nm, set = rownames(r), r)
}))
cam[, genotype := sub("_.*", "", contrast)][, contrast := sub("^[0-9]+_", "", contrast)]
cam <- merge(cam, sm[, .(genotype = as.character(genotype), contrast, status)], by = c("genotype", "contrast"), all.x = TRUE)
fwrite(cam[order(genotype, contrast, PValue)], file.path(OUT, "camera_literature_sets.tsv"), sep = "\t")

# our memory classes vs literature classes (tillering)
ov <- list()
for (g in c("5000", "8008")) {
  f <- file.path(DE, sprintf("memory_classes_%s_tillering.tsv.gz", g)); if (!file.exists(f)) next
  mc <- fread(f)
  for (sp in c("Athaliana", "Zmays")) for (cl in c("++", "--", "+-", "-+")) {
    litg <- unique(map[lit_species == sp & class == cl, group_id])
    for (ours in setdiff(unique(mc$class), "none")) {
      a <- mc$class == ours; b <- mc$id %in% litg
      ft <- fisher.test(table(factor(a, c(TRUE, FALSE)), factor(b, c(TRUE, FALSE))), alternative = "greater")
      ov[[length(ov) + 1]] <- data.table(genotype = g, our_class = ours, lit_species = sp, lit_class = cl,
        n_ours = sum(a), n_lit_tested = sum(b), overlap = sum(a & b), odds_ratio = unname(ft$estimate), p = ft$p.value,
        groups = paste(head(mc$id[a & b], 30), collapse = ","))
    }
  }
}
ov <- rbindlist(ov); ov[, fdr := p.adjust(p, "BH"), by = genotype]
fwrite(ov[order(genotype, p)], file.path(OUT, "memory_class_overlap_literature.tsv"), sep = "\t")

# heatmap of curated candidates
cg <- unique(map[grepl("^curated", class), .(group_id, label = paste0(symbol, " | ", sub("^SB:", "", group_id)), class)])
for (g in c("5000", "8008")) {
  fs <- files[grepl(sprintf("de_%s_", g), basename(files))]
  m <- sapply(fs, function(f) { dt <- fread(f); dt$log2FC_shrunk[match(cg$group_id, dt$id)] })
  colnames(m) <- sub(sprintf("^de_%s_", g), "", sub("\\.tsv\\.gz$", "", basename(fs))); rownames(m) <- make.unique(cg$label)
  m <- m[rowSums(!is.na(m)) > 0, , drop = FALSE]; m[is.na(m)] <- 0
  fwrite(data.table(group = rownames(m), m), file.path(OUT, sprintf("curated_lfc_%s.tsv", g)), sep = "\t")
  ann <- data.frame(category = sub("curated:", "", cg$class[match(rownames(m), make.unique(cg$label))]), row.names = rownames(m))
  br <- seq(-4, 4, length.out = 101)
  pheatmap(pmax(pmin(m, 4), -4), breaks = br, cluster_cols = FALSE, annotation_row = ann, fontsize_row = 5,
           main = paste("Curated candidates, shrunk log2FC,", g), filename = file.path(OUT, sprintf("curated_heatmap_%s.png", g)),
           width = 10, height = max(6, nrow(m) * 0.09))
}
cat("done\n")
