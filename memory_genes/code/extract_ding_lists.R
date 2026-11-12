# Extract memory-class gene lists from Ding et al. 2013 (Arabidopsis, BMC Plant Biol 13:229, Additional file 4)
# and Ding et al. 2014 (maize, BMC Plant Biol 14:141, Additional file 4).
# Classes (sheet names): ++ -- +- -+ (memory) and += -= (non-memory). Each sheet has one row per gene with
# water / S1 / S3 FPKM and log2FC; the maize sheets also start with a small summary table of Arabidopsis homologs.
suppressMessages(library(readxl))
D <- "/Storage/data1/jorge.munoz/DOLORES/memory_genes/literature"
take <- function(f, species, idre, symcol = TRUE) {
  out <- list()
  for (s in excel_sheets(f)) {
    if (!s %in% c("++", "--", "+-", "-+", "+=", "-=")) next
    d <- suppressMessages(read_excel(f, sheet = s, col_names = FALSE, .name_repair = "minimal"))
    d <- as.data.frame(d)
    hdr <- which(apply(d, 1, function(r) any(r %in% "gene_id")))[1]
    if (is.na(hdr)) next
    cn <- as.character(unlist(d[hdr, ])); x <- d[(hdr + 1):nrow(d), , drop = FALSE]; names(x) <- cn
    gcol <- if (species == "Athaliana") "gene_id" else "gene"
    ids <- as.character(x[[gcol]]); keep <- grepl(idre, ids)
    num <- function(n) suppressWarnings(as.numeric(x[[n]][keep]))
    out[[s]] <- data.frame(species = species, class = s, gene_id = ids[keep],
                           symbol = if (species == "Athaliana") as.character(x[["gene"]][keep]) else NA,
                           fpkm_water = num("water"), fpkm_S1 = num("S1"), fpkm_S3 = num("S3"),
                           lfc_S1_water = num("S1_water log2(fold_change)"), lfc_S3_S1 = num("S3_S1 log2(fold_change)"))
  }
  do.call(rbind, out)
}
a <- take(file.path(D, "PMC3879431/1471-2229-13-229-S4.xls"), "Athaliana", "^AT[1-5CM]G[0-9]+")
z <- take(file.path(D, "PMC4081654/1471-2229-14-141-S4.xlsx"), "Zmays_v3", "^(GRMZM|AC[0-9]+\\.[0-9]+_FG|EF|AY)")
r <- rbind(a, z)
r$source <- ifelse(r$species == "Athaliana", "Ding et al. 2013 BMC Plant Biol 13:229", "Ding et al. 2014 BMC Plant Biol 14:141")
write.table(r, "/Storage/data1/jorge.munoz/DOLORES/memory_genes/ding_memory_classes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(table(r$species, r$class))
