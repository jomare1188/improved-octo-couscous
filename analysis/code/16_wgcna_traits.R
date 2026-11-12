#!/usr/bin/env Rscript
# Which leaf co-expression modules follow drought HISTORY rather than current drought?
# Each module eigengene (ME, one value per library) is modelled with limma:
#   ME ~ genotype + propagule + maturation + late_run + current + repeated + maternal
#     current  = under drought when sampled (D1, D3, propagule D1)
#     repeated = 3rd drought (Mt_D3, Mm_D3) on top of the current-drought effect: the plant's own history
#     maternal = propagule of a stressed mother (Ft, Fm), on top of current drought: the mother's history
#     late_run = Mt_W, the only group sequenced on runs 61/63/73 (absorbs that run and the age difference)
# A "history module" has FDR < 0.05 for repeated or maternal. The same model is fitted in each genotype to
# check that the effect has the same sign in both. Then: GO terms (topGO), overlap with memory classes and
# DEG lists, hub groups, and module scores in buds (T vs M is a clean comparison, W is run-confounded).
# Run from analysis/:  qsub -N wgtraits code/run_r.sh code/16_wgcna_traits.R
source("code/common.R")
suppressMessages({ library(limma); library(topGO) })
OUT <- file.path(AN, "results/wgcna")
P <- plot_setup()
mods <- fread(file.path(OUT, "modules.tsv")); MEs <- fread(file.path(OUT, "eigengenes.tsv"))
tr <- fread(file.path(OUT, "traits.tsv"), colClasses = list(character = "genotype"))
tr <- tr[match(MEs$sample, sample)]
me <- t(as.matrix(MEs[, -1])); colnames(me) <- MEs$sample
terms <- c("current", "repeated", "maternal", "late_run", "propagule", "maturation")

fit_model <- function(idx, with_genotype) {
  f <- if (with_genotype) ~ g8008 + propagule + maturation + late_run + current + repeated + maternal
       else ~ propagule + maturation + late_run + current + repeated + maternal
  X <- model.matrix(f, tr[idx]); X <- X[, colSums(abs(X)) > 0, drop = FALSE]
  fit <- eBayes(lmFit(me[, idx], X))
  rbindlist(lapply(intersect(terms, colnames(X)), function(k) {
    tt <- topTable(fit, coef = k, number = Inf, sort.by = "none")
    data.table(module = rownames(tt), term = k, effect = tt$logFC, t = tt$t, P = tt$P.Value, FDR = tt$adj.P.Val)
  }))
}
pooled <- fit_model(seq_len(nrow(tr)), TRUE)[, model := "both genotypes"]
per_g <- rbindlist(lapply(c("8008", "5000"), function(g) fit_model(which(tr$genotype == g), FALSE)[, model := g]))
mt <- rbind(pooled, per_g)
size <- mods[, .N, by = module]
mt <- merge(mt, size, by = "module")
fwrite(mt, file.path(OUT, "module_trait_tests.tsv"), sep = "\t")

hist_mods <- pooled[term %in% c("repeated", "maternal") & FDR < 0.05 & module != "M00", unique(module)]
# modules with a history signal in at least one genotype (weaker evidence; M00 = unassigned, left out)
any_mods <- mt[term %in% c("repeated", "maternal") & FDR < 0.05 & module != "M00", unique(module)]
cons <- dcast(mt[module %in% any_mods & term %in% c("current", "repeated", "maternal")], module + term ~ model, value.var = "t")
fwrite(cons, file.path(OUT, "history_modules_consistency.tsv"), sep = "\t")
cat("history modules (pooled):", hist_mods, " in one genotype:", any_mods, "\n"); print(cons)

# ---------------------------------------------------------------- module contents
lab <- group_labels()
kme <- fread(file.path(OUT, "kme_all.tsv.gz"))
hubs <- rbindlist(lapply(setdiff(unique(mods$module), "M00"), function(m) {
  x <- merge(mods[module == m], lab[, .(group_id, name, description)], by = "group_id", all.x = TRUE)
  x[order(-kME_own)][1:min(15, .N)][, rank := .I]
}))
fwrite(hubs, file.path(OUT, "module_hubs.tsv"), sep = "\t")

universe <- mods$group_id
g2go <- readMappings(file.path(ORTHO, "group2go_list.tsv"))
uni_go <- intersect(universe, names(g2go))
go <- rbindlist(lapply(setdiff(unique(mods$module), "M00"), function(m) {
  gl <- factor(as.integer(uni_go %in% mods[module == m, group_id]), levels = 0:1); names(gl) <- uni_go
  if (sum(gl == 1) < 5) return(NULL)
  tg <- suppressMessages(new("topGOdata", ontology = "BP", allGenes = gl, annot = annFUN.gene2GO, gene2GO = g2go, nodeSize = 10))
  rs <- suppressMessages(runTest(tg, algorithm = "weight01", statistic = "fisher"))
  tb <- as.data.table(GenTable(tg, p = rs, topNodes = 15, numChar = 120))
  tb[, p := as.numeric(sub("<", "", p))][p < 0.01][, module := m]
}))
fwrite(go, file.path(OUT, "module_go_bp.tsv"), sep = "\t")

# overlap with memory classes and DEG lists (Fisher, universe = network groups)
lists <- list()
for (g in c("8008", "5000")) {
  mc <- fread(file.path(AN, sprintf("results/de/group/memory_classes_%s_tillering.tsv.gz", g)))
  for (cl in c("++", "--", "0+", "0-")) lists[[sprintf("%s memory %s", g, cl)]] <- mc[class == cl, id]
  for (k in c("Mt_D3_vs_Mt_D1", "Ft_D1_vs_Fw_D1", "Mt_D1_vs_Mt_W")) {
    x <- rd_de(g, k); lists[[sprintf("%s %s up", g, k)]] <- x[deg == "up", id]; lists[[sprintf("%s %s down", g, k)]] <- x[deg == "down", id]
  }
}
ov <- rbindlist(lapply(setdiff(unique(mods$module), "M00"), function(m) rbindlist(lapply(names(lists), function(l) {
  a <- universe %in% mods[module == m, group_id]; b <- universe %in% lists[[l]]
  if (sum(b) < 5) return(NULL)
  ft <- fisher.test(table(factor(a, c(TRUE, FALSE)), factor(b, c(TRUE, FALSE))), alternative = "greater")
  data.table(module = m, list = l, module_size = sum(a), list_in_network = sum(b), overlap = sum(a & b),
             odds_ratio = round(unname(ft$estimate), 2), p = ft$p.value)
}))))
ov[, FDR := p.adjust(p, "BH")]
fwrite(ov[order(p)], file.path(OUT, "module_list_overlaps.tsv"), sep = "\t")

# ---------------------------------------------------------------- module scores in buds
xb <- readRDS(file.path(OUT, "bud_expr_centred.rds"))
bmeta <- read_meta(); bmeta <- bmeta[colnames(xb), ]
z <- t(scale(t(xb)))
bs <- rbindlist(lapply(setdiff(unique(mods$module), "M00"), function(m) {
  gs <- intersect(mods[module == m, group_id], rownames(z))
  data.table(module = m, sample = colnames(z), genotype = bmeta$genotype, condition = bmeta$condition, score = colMeans(z[gs, , drop = FALSE]))
}))
fwrite(bs, file.path(OUT, "bud_module_scores.tsv"), sep = "\t")
btest <- bs[condition %in% c("T", "M"), {
  a <- score[condition == "T"]; b <- score[condition == "M"]
  .(T_minus_M = round(mean(a) - mean(b), 3), p_T_vs_M = if (length(a) > 1 && length(b) > 1) t.test(a, b)$p.value else NA_real_)
}, by = .(module, genotype)]
fwrite(btest, file.path(OUT, "bud_module_T_vs_M.tsv"), sep = "\t")

# ---------------------------------------------------------------- figures
pm <- mt[module != "M00"][, model := factor(model, c("both genotypes", "8008", "5000"))][, term := factor(term, terms, c("current drought", "repeated (D3)", "maternal stress", "late run (Mt_W)", "propagule", "maturation"))]
ord <- pm[term == "current drought" & model == "both genotypes"][order(t), module]
pm[, module := factor(module, ord)]
pm[, star := fifelse(FDR < 0.01, "**", fifelse(FDR < 0.05, "*", ""))]
p1 <- ggplot(pm, aes(term, module, fill = pmax(pmin(t, 10), -10))) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = star, colour = abs(t) > 6), size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = P$ink)) +
  facet_wrap(~ model, nrow = 1) +
  scale_fill_gradientn(colours = P$div, limits = c(-10, 10), name = "t (capped ±10)") +
  scale_y_discrete(labels = function(m) sprintf("%s (%s)", m, size$N[match(m, size$module)])) +
  labs(x = NULL, y = "Module (number of groups)", title = "What does each leaf module follow?",
       subtitle = "limma on module eigengenes, both genotypes together and each alone. * FDR < 0.05, ** FDR < 0.01.\nRed = higher with the factor, blue = lower. Unassigned groups (M00) not shown.") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p1, file.path(OUT, "module_trait.png"), 14, max(5, 0.3 * length(ord) + 2.5))

plot_mods <- if (length(hist_mods)) hist_mods else any_mods
if (length(plot_mods)) {
  ev <- melt(MEs[, c("sample", plot_mods), with = FALSE], id.vars = "sample", variable.name = "module", value.name = "ME")
  ev <- merge(ev, tr[, .(sample, genotype, condition)], by = "sample")
  lv <- c("Mt_W", "Mt_D1", "Mt_D3", "Mm_W", "Mm_D1", "Mm_D3", "Fw_W", "Fw_D1", "Ft_D1", "Fm_D1")
  ev[, condition := factor(condition, lv)]
  p2 <- ggplot(ev, aes(condition, ME, colour = genotype)) +
    geom_hline(yintercept = 0, colour = P$grid) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    facet_wrap(~ module, scales = "free_y") + scale_colour_manual(values = P$geno, name = "Genotype") +
    labs(x = NULL, y = "Module eigengene", title = "Modules that follow drought history",
         subtitle = "One point per library. Mother leaves: Mt (tillering), Mm (maturation). Propagules: Fw, Ft, Fm.") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  if (!length(hist_mods)) p2 <- p2 + labs(title = "Modules with a history signal in one genotype only")
  save_plot(p2, file.path(OUT, "history_modules.png"), max(6, 3.6 * min(3, length(plot_mods))), 3 + 2.6 * ceiling(length(plot_mods) / 3))
  bsh <- bs[module %in% plot_mods][, condition := factor(condition, c("W", "T", "M"))]
  p3 <- ggplot(bsh, aes(condition, score, colour = condition)) + geom_point(size = 2.2) +
    facet_grid(genotype ~ module) +
    scale_colour_manual(values = setNames(unname(P$hist), c("W", "T", "M")), guide = "none") +
    labs(x = "Bud group (W on another run)", y = "Module score (mean z)", title = "History modules in buds")
  save_plot(p3, file.path(OUT, "history_modules_in_buds.png"), max(6, 1.6 * length(plot_mods) + 2), 4.5)
}
cat("done\n")
