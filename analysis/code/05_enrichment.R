#!/usr/bin/env Rscript
# GO (topGO, weight01 + Fisher, BP and MF) and KEGG pathway (clusterProfiler::enricher) enrichment at the
# allele-group level, for up and down DEGs of every contrast and for the tillering memory classes.
# Universe = groups tested in that contrast (non-NA padj). Group GO = union over member alleles;
# group KEGG pathways = union of eggNOG pathways of member alleles.
suppressMessages({ library(data.table); library(topGO); library(clusterProfiler) })
D  <- "/Storage/data1/jorge.munoz/DOLORES"
DE <- file.path(D, "analysis/results/de/group")
OUT <- file.path(D, "analysis/results/enrichment"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

g2go <- readMappings(file.path(D, "orthology/allele_groups/group2go_list.tsv"))
gene2group <- fread(file.path(D, "orthology/allele_groups/gene2group.tsv"))
ann <- fread(file.path(D, "annotation/go/gene_annotation.tsv"), select = c("gene_id", "kegg_pathway"))[kegg_pathway != ""]
kp <- merge(ann, gene2group, by = "gene_id")[, .(pathway = unlist(strsplit(kegg_pathway, ","))), by = group_id]
kp <- unique(kp[, .(pathway, group_id)])
kn <- fread(file.path(D, "analysis/data/kegg_pathway_names.tsv"), header = FALSE, col.names = c("pathway", "name"))

run_go <- function(sel, universe, tag) {
  universe <- intersect(universe, names(g2go)); if (sum(universe %in% sel) < 5) return(NULL)
  rbindlist(lapply(c("BP", "MF"), function(ont) {
    gl <- factor(as.integer(universe %in% sel), levels = 0:1); names(gl) <- universe
    tg <- suppressMessages(new("topGOdata", ontology = ont, allGenes = gl, annot = annFUN.gene2GO, gene2GO = g2go, nodeSize = 10))
    rs <- suppressMessages(runTest(tg, algorithm = "weight01", statistic = "fisher"))
    tb <- GenTable(tg, p = rs, topNodes = min(200, length(score(rs))), numChar = 200)
    tb$p <- as.numeric(sub("<", "", tb$p)); tb <- tb[tb$p < 0.01, ]
    if (!nrow(tb)) return(NULL)
    data.table(set = tag, ontology = ont, tb)
  }))
}
run_kegg <- function(sel, universe, tag) {
  e <- enricher(intersect(sel, universe), universe = universe, TERM2GENE = kp, TERM2NAME = kn,
                pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10)
  if (is.null(e) || !nrow(as.data.frame(e))) return(NULL)
  data.table(set = tag, as.data.frame(e))
}

go_all <- list(); kg_all <- list()
for (f in list.files(DE, "^de_.*\\.tsv\\.gz$", full.names = TRUE)) {
  dt <- fread(f); nm <- sub("^de_", "", sub("\\.tsv\\.gz$", "", basename(f)))
  uni <- dt[!is.na(padj), id]
  for (dir in c("up", "down")) {
    sel <- dt[deg == dir, id]; tag <- paste(nm, dir, sep = ":")
    cat(tag, length(sel), "\n")
    go_all[[tag]] <- run_go(sel, uni, tag); kg_all[[tag]] <- run_kegg(sel, uni, tag)
  }
}
for (g in c("5000", "8008")) {
  f <- file.path(DE, sprintf("memory_classes_%s_tillering.tsv.gz", g)); if (!file.exists(f)) next
  mc <- fread(f)[!is.na(D1_vs_W_padj)]
  for (cl in setdiff(unique(mc$class), c("none", "+=", "-="))) {
    tag <- sprintf("%s_memory:%s", g, cl); sel <- mc[class == cl, id]
    cat(tag, length(sel), "\n")
    go_all[[tag]] <- run_go(sel, mc$id, tag); kg_all[[tag]] <- run_kegg(sel, mc$id, tag)
  }
}
fwrite(rbindlist(go_all, fill = TRUE), file.path(OUT, "go_enrichment.tsv"), sep = "\t")
fwrite(rbindlist(kg_all, fill = TRUE), file.path(OUT, "kegg_enrichment.tsv"), sep = "\t")
cat("done\n")
