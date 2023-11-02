library(rlang, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(graph, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(GO.db, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(SparseM, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(topGO, lib.loc = "/Storage/data1/jorge.munoz/TOPGO/data")
library(dplyr, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(Rgraphviz, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(ggplot2, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(scales, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(clusterProfiler, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")

setwd("/Storage/data1/jorge.munoz/DOLORES/functional_enrichmnent/go_enrichment/code")

GOenrichment <- function(contrast, ontology) {

title <- switch(ontology,
  "BP" = "GO Biological process",
  "CC" = "GO Cellular component",
  "MF" = "GO Molecular function",
  "Valor de ontology no reconocido"
)

geneID2GO <- readMappings(file = "../data/all_GO_for_enrichment.tsv")

class_table <- read.table("../data/all_diff_expr_classification.tsv", header = T, sep = ",")
geneNames <- read.table("../data/diff_expr.ids", header = F)

myInterestingGenes <-class_table[class_table$Group==contrast,1]
geneList <- as.numeric(as.integer(geneNames$V1 %in% myInterestingGenes))
names(geneList) <- geneNames$V1
geneList <- as.factor(geneList)

GOdata <- new("topGOdata",
              ontology = ontology,
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

allGO <- usedGO(GOdata)
Classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic')
# Filter not significant values for classic algorithm
table1 <- filter(table, Classic < 0.05 )
# Performing BH correction on our p values FDR
p.adj <- round(p.adjust(table1$Classic,method="BH"),digits = 4)
# Create the file with all the statistics from GO analysis
all_res_final <<- cbind(table1,p.adj)
all_res_final <<- all_res_final[order(all_res_final$p.adj),]
# Get list of significant GO before multiple testing correction
results.table.p = all_res_final[which(all_res_final$Classic <=0.05),]
# Get list of significant GO after multiple testing correction
results.table.bh = all_res_final[which(all_res_final$p.adj<=0.05),]
# Save first top 50 ontolgies sorted by adjusted pvalues
write.table(results.table.bh, file = paste0("../results/", contrast, "_", ontology, ".csv"), quote=FALSE, row.names=FALSE, sep = ",")
ntop <- 12
ggdata <- all_res_final[1:ntop,]

ggdata <- ggdata[complete.cases(ggdata), ]

aux <- go2term(all_res_final$GO.ID)
colnames(aux) <- c("GO.ID", "Lterm")

ggdata <- merge(ggdata, aux, by = "GO.ID")

ggdata$Classic <- as.numeric(ggdata$Classic)

ggdata <- ggdata[order(ggdata$Classic),]
ggdata$Lterm <- factor(ggdata$Lterm, levels = rev(ggdata$Lterm)) # fixes order

gg1 <- ggplot(ggdata, aes(x = Lterm, y = -log10(Classic)))+
  geom_point(size = 6, colour = "black") +
  scale_size(range = c(2.5,12.5)) +
  xlab('GO Term') +
  ylab('-log10(p)') +
  labs(title = title) +
  coord_flip() +
  theme_bw(base_size = 24)

ggsave(paste0("../results/", contrast, "_", ontology, ".png"), device = "png", width = 60, height = 30, dpi = 300, units = "cm")
}
# 5_M_D_vs_5_M_W
GOenrichment("5_M_D_vs_5_M_W_up", "BP")
GOenrichment("5_M_D_vs_5_M_W_down", "BP")

GOenrichment("5_M_D_vs_5_M_W_up", "CC")
GOenrichment("5_M_D_vs_5_M_W_down", "CC")

GOenrichment("5_M_D_vs_5_M_W_up", "MF")
GOenrichment("5_M_D_vs_5_M_W_down", "MF")

# 5_M_D_vs_8_M_D
GOenrichment("5_M_D_vs_8_M_D_up", "BP")
GOenrichment("5_M_D_vs_8_M_D_down", "BP")

GOenrichment("5_M_D_vs_8_M_D_up", "CC")
GOenrichment("5_M_D_vs_8_M_D_down", "CC")

GOenrichment("5_M_D_vs_8_M_D_up", "MF")
GOenrichment("5_M_D_vs_8_M_D_down", "MF")
# 5_M_W_vs_8_M_W
GOenrichment("5_M_W_vs_8_M_W_up", "BP")
GOenrichment("5_M_W_vs_8_M_W_down", "BP")

GOenrichment("5_M_W_vs_8_M_W_up", "CC")
GOenrichment("5_M_W_vs_8_M_W_down", "CC")

GOenrichment("5_M_W_vs_8_M_W_up", "MF")
GOenrichment("5_M_W_vs_8_M_W_down", "MF")
# 5_T_D_vs_8_T_D
GOenrichment("5_T_D_vs_8_T_D_up", "BP")
GOenrichment("5_T_D_vs_8_T_D_down", "BP")

GOenrichment("5_T_D_vs_8_T_D_up", "CC")
GOenrichment("5_T_D_vs_8_T_D_down", "CC")

GOenrichment("5_T_D_vs_8_T_D_up", "MF")
GOenrichment("5_T_D_vs_8_T_D_down", "MF")
# 8_M_D_vs_8_M_W
GOenrichment("8_M_D_vs_8_M_W_up", "BP")
GOenrichment("8_M_D_vs_8_M_W_down", "BP")

GOenrichment("8_M_D_vs_8_M_W_up", "CC")
GOenrichment("8_M_D_vs_8_M_W_down", "CC")

GOenrichment("8_M_D_vs_8_M_W_up", "MF")
GOenrichment("8_M_D_vs_8_M_W_down", "MF")

