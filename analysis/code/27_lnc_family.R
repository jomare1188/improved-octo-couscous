# Genome-wide lncRNA search, last step: the lncRNAs that differ between T and M buds in both genotypes,
# the sequence family most of them belong to, and the figure for report_story.md (section 7).
#   cd analysis && bash code/run_r.sh code/27_lnc_family.R      (after 26_lnc_de.R)
# Family: shared T-vs-M lncRNAs are compared all against all (blastn, >= 150 nt); the largest connected set is
# the core; every lncRNA transcript matching a core member (>= 200 nt, >= 80% identity) joins the family.
# Identity of the family (NCBI core_nt, nr, RepBase, Rfam, pan-transcriptome) was searched by hand: see
# results/lncrna_genome/family_identity.md.
source("code/common.R"); suppressMessages(library(patchwork)); st <- plot_setup()
W   <- "/Storage/data2/jorge.munoz/lnc_genome"
OUT <- file.path(AN, "results/lncrna_genome")
BLASTN <- "/Storage/progs/miniconda3/bin/blastn"
tmp <- file.path(W, "family"); dir.create(tmp, showWarnings = FALSE)

de  <- fread(file.path(OUT, "lnc_de_all_contrasts.tsv.gz"), colClasses = list(character = "genotype"))
cl  <- fread(file.path(OUT, "lncRNA_clusters.tsv"))
t2g <- fread(file.path(OUT, "lnc_tx2group.tsv"))
ca  <- fread(file.path(OUT, "candidates.tsv"))

# ---- lncRNAs different between T and M buds in both genotypes, same direction
tm <- merge(de[genotype == "8008" & contrast == "T_vs_M"], de[genotype == "5000" & contrast == "T_vs_M"], by = "id", suffixes = c("_8008", "_5000"))
shared <- tm[deg_8008 != "" & deg_8008 == deg_5000]
expected <- sum(tm$deg_8008 != "") * sum(tm$deg_5000 != "") / nrow(tm)
cat(sprintf("T vs M lncRNAs: 8008 %d, 5000 %d, tested in both %d, shared same direction %d (expected by chance %.2f)\n",
            sum(tm$deg_8008 != ""), sum(tm$deg_5000 != ""), nrow(tm), nrow(shared), expected))

# ---- sequence family
fa <- readLines(file.path(W, "lnc.fa")); h <- grep("^>", fa)
seqs <- data.table(tx = sub("^>(\\S+).*", "\\1", fa[h]), from = h + 1, to = c(h[-1] - 1, length(fa)))
write_fa <- function(ids, f) { s <- seqs[tx %in% ids]
  writeLines(unlist(lapply(seq_len(nrow(s)), function(i) c(paste0(">", s$tx[i]), toupper(fa[s$from[i]:s$to[i]])))), f) }
q <- file.path(tmp, "shared.fa"); write_fa(sub("^LNC:", "", shared$id), q)
bl <- function(query, subject, minlen, minid = 0) {
  x <- fread(cmd = sprintf("%s -query %s -subject %s -evalue 1e-10 -max_target_seqs 5000 -outfmt '6 qseqid sseqid pident length' 2>/dev/null", BLASTN, query, subject), header = FALSE)
  x[V4 >= minlen & V3 >= minid]
}
aa <- bl(q, q, 150)[V1 != V2]
comp <- setNames(seq_len(nrow(shared)), sub("^LNC:", "", shared$id))
repeat { old <- comp; for (i in seq_len(nrow(aa))) { m <- min(comp[aa$V1[i]], comp[aa$V2[i]]); comp[aa$V1[i]] <- m; comp[aa$V2[i]] <- m }
  if (identical(old, comp)) break }
core <- names(comp)[comp == as.integer(names(which.max(table(comp))))]
qc <- file.path(tmp, "core.fa"); write_fa(core, qc)
fam_tx <- unique(bl(qc, file.path(W, "lnc.fa"), 200, 80)$V2)
fam <- unique(t2g[transcript_id %in% fam_tx, group_id])
cat(sprintf("core %d of %d shared; family: %d transcripts, %d clusters, %d genomic copies\n",
            length(core), nrow(shared), length(fam_tx), length(fam), sum(cl[group_id %in% fam, n_copies])))

wide <- dcast(de[id %in% fam & contrast %in% c("T_vs_M", "Mt_D1_vs_Mt_W", "Mt_D3_vs_Mt_W")], id ~ genotype + contrast, value.var = c("log2FC", "padj"))
famtab <- merge(cl[, .(id = group_id, n_copies, max_exons, length_representative, max_tpm)], wide, by = "id")
famtab[, core := sub("^LNC:", "", id) %in% core]
fwrite(famtab[order(-max_tpm)], file.path(OUT, "bud_family_clusters.tsv"), sep = "\t")
sh <- merge(shared[, .(id, log2FC_8008, padj_8008, log2FC_5000, padj_5000, direction = fifelse(deg_8008 == "up", "higher in T", "higher in M"))],
            cl[, .(id = group_id, n_copies, max_exons, length_representative, max_tpm)], by = "id")
sh[, family := id %in% fam]
fwrite(sh[order(-family, id)], file.path(OUT, "bud_shared_lncRNAs.tsv"), sep = "\t")
print(de[id %in% fam & contrast %in% c("T_vs_M", "Mt_D1_vs_Mt_W", "Mt_D3_vs_Mt_W"),
         .(tested = sum(!is.na(padj)), median_log2FC = round(median(log2FC[!is.na(padj)]), 2), up = sum(deg == "up"), down = sum(deg == "down")), by = .(genotype, contrast)])

# ---- family expression per library (sum of TPM over the family's transcripts)
meta <- as.data.table(read_meta())
tpm <- rbindlist(lapply(seq_len(nrow(meta)), function(i) {
  x <- fread(file.path(W, "quant", meta$library[i], "quant.sf"), select = c("Name", "TPM"))
  data.table(sample = meta$sample[i], family_tpm = sum(x[Name %in% fam_tx, TPM]))
}))
tpm <- merge(tpm, meta[, .(sample, genotype, set, condition)], by = "sample")
fwrite(tpm[order(set, genotype, condition)], file.path(OUT, "family_tpm_by_library.tsv"), sep = "\t")

# ---- figure: the filters (left) and the family in buds (right)
steps <- data.table(
  label = c("Transcripts assembled\nfrom the 48 libraries",
            "New, outside annotated genes,\n200 nt or longer",
            "No coding signal: no ORF of 100 aa,\nno protein match, CPC2 non-coding,\nnot beside a gene on the same strand",
            "Clusters of copies\n(90% identity or more)",
            "High confidence: 2 or more exons,\nexpressed in 2 or more libraries"),
  n = c(as.integer(sub(" .*", "", grep("consensus transcripts written", readLines(file.path(OUT, "cmp.stats")), value = TRUE))),
        nrow(ca), sum(ca$lncRNA == 1), nrow(cl), sum(cl$high_confidence)))
steps[, label := factor(label, levels = rev(label))]
pA <- ggplot(steps, aes(n, label)) +
  geom_col(fill = "#B8B8B8", width = 0.62) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.12, size = 3.3, colour = st$ink) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "A. From reads to lncRNA candidates", x = "Transcripts or clusters", y = NULL) +
  theme(axis.text.y = element_text(colour = st$ink, size = 8.2, lineheight = 0.95), axis.text.x = element_blank(),
        panel.grid.major = element_blank(), plot.title = element_text(size = 10.5))

pb <- tpm[set == "bud"]
lev <- c(W = "W\nwatered", T = "T\nstressed\n135 d ago", M = "M\nstressed\n6 d ago")
pb[, step := factor(lev[condition], levels = lev)][, y := log2(family_tpm + 1)][, genotype := factor(genotype, levels = c("8008", "5000"))]
pm <- pb[, .(y = mean(y)), by = .(genotype, step)]
cols <- setNames(unname(st$hist), lev)
pB <- ggplot(pb, aes(step, y)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#F2F2F2") +
  annotate("text", x = 1, y = Inf, label = "other run", vjust = 1.6, size = 2.8, colour = st$muted) +
  geom_errorbar(data = pm, aes(ymin = y, ymax = y, colour = step), width = 0.45, linewidth = 1.1) +
  geom_point(aes(fill = step), shape = 21, colour = "white", size = 3.2, stroke = 0.6, position = position_jitter(width = 0.08, height = 0, seed = 4)) +
  scale_x_discrete() + scale_fill_manual(values = cols, guide = "none") + scale_colour_manual(values = cols, guide = "none") +
  facet_wrap(~ genotype) +
  labs(title = sprintf("B. One family of %d clusters is lower in T than in M buds", length(fam)),
       x = NULL, y = "log2 (summed TPM of the family + 1)") +
  theme(axis.text.x = element_text(colour = st$ink, size = 8.2), panel.grid.major.x = element_blank(), plot.title = element_text(size = 10.5))
p <- pA + pB + plot_layout(widths = c(1, 1.15))
save_plot(p, file.path(AN, "results/plots/story/lncrna.png"), 11, 4.6)
cat("done\n")
