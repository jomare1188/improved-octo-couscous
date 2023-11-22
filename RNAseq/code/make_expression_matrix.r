library("dplyr", lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library("dynamicTreeCut", lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library("fastcluster", lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library("WGCNA", lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library("igraph", lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")

setwd("/Storage/data1/jorge.munoz/DOLORES/RNAseq/code")

## make histogram for DEG in all contraste (upregulated and downregulated)
contrast_results_dir <- "/Storage/data1/jorge.munoz/DOLORES/RNAseq/results"
log2foldthreshold <- 2
#contrast 1
library("dplyr", lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")

#old
#contrast_1 <-as.data.frame(read.table(paste0(contrast_results_dir,"/","5_M_D_vs_8_M_D", "/", "5_M_D_vs_8_M_D.DESeq.DE_results"), header = T) %>% filter(abs(log2FoldChange) > log2foldthreshold ) %>% filter((padj) < 0.01),header = F)

#contrast_2 <-as.data.frame(read.table(paste0(contrast_results_dir,"/","5_T_D_vs_8_T_D", "/", "5_T_D_vs_8_T_D.DESeq.DE_results"), header = T) %>% filter(abs(log2FoldChange) > log2foldthreshold ) %>% filter((padj) < 0.01),header = F)

#contrast_3 <-as.data.frame(read.table(paste0(contrast_results_dir,"/","5_M_W_vs_8_M_W", "/", "5_M_W_vs_8_M_W.DESeq.DE_results"), header = T) %>% filter(abs(log2FoldChange) > log2foldthreshold ) %>% filter((padj) < 0.01),header = F)

#contrast_4 <-as.data.frame(read.table(paste0(contrast_results_dir,"/","5_M_D_vs_5_M_W", "/", "5_M_D_vs_5_M_W.DESeq.DE_results"), header = T) %>% filter(abs(log2FoldChange) > log2foldthreshold ) %>% filter((padj) < 0.01),header = F)

#contrast_5 <-as.data.frame(read.table(paste0(contrast_results_dir,"/","8_M_D_vs_8_M_W", "/", "8_M_D_vs_8_M_W.DESeq.DE_results"), header = T) %>% filter(abs(log2FoldChange) > log2foldthreshold ) %>% filter((padj) < 0.01),header = F)

# new

contrast_1 <-as.data.frame(read.table(paste0(contrast_results_dir,"/","5_M_W_vs_5_T_D", "/", "5_M_W_vs_5_T_D.DESeq.DE_results"), header = T) %>% filter(abs(log2FoldChange) > log2foldthreshold ) %>% filter((padj) < 0.01),header = F)

contrast_2 <-as.data.frame(read.table(paste0(contrast_results_dir,"/","8_M_W_vs_8_T_D", "/", "8_M_W_vs_8_T_D.DESeq.DE_results"), header = T) %>% filter(abs(log2FoldChange) > log2foldthreshold ) %>% filter((padj) < 0.01),header = F)




df3 <- data.frame(matrix(ncol = 3, nrow = 4))
colnames(df3) <- c("Contrast", "Type", "DEGs") 
# contrastast 5_M_D_vs_5_T_D
# up
df3[1,1] <- "5_M_W_vs_5_T_D"
df3[1,2] <- "Up" 
df3[1,3] <- dim(contrast_1[contrast_1$log2FoldChange > 1,])[1]
names <- rownames(contrast_1[contrast_1$log2FoldChange > 1,])
write.table(names, paste0(contrast_results_dir, "/","5_M_W_vs_5_T_D", "/", "5_M_W_vs_5_T_D_up"), col.names = F, row.names = F, quote = F)
# down
df3[2,1] <- "5_M_W_vs_5_T_D"
df3[2,2] <- "Down" 
df3[2,3] <- dim(contrast_1[contrast_1$log2FoldChange < 1,])[1]
names <- rownames(contrast_1[contrast_1$log2FoldChange < 1,])
write.table(names, paste0(contrast_results_dir, "/","5_M_W_vs_5_T_D", "/","5_M_W_vs_5_T_D_down"), col.names = F, row.names = F, quote = F)

# contrastast 8_M_D_vs_8_T_D
# up
df3[3,1] <- "8_M_W_vs_8_T_D"
df3[3,2] <- "Up"
df3[3,3] <- dim(contrast_2[contrast_2$log2FoldChange > 1,])[1]
names <- rownames(contrast_2[contrast_2$log2FoldChange > 1,])
write.table(names, paste0(contrast_results_dir, "/","8_M_W_vs_8_T_D", "/", "8_M_W_vs_8_T_D_up"), col.names = F, row.names = F, quote = F)
# down 
df3[4,1] <- "8_M_W_vs_8_T_D"
df3[4,2] <- "Down"
df3[4,3] <- dim(contrast_2[contrast_2$log2FoldChange < 1,])[1]
names <- rownames(contrast_2[contrast_2$log2FoldChange < 1,])
write.table(names, paste0(contrast_results_dir, "/","8_M_W_vs_8_T_D", "/","8_M_W_vs_8_T_D_down"), col.names = F, row.names = F, quote = F)

# contrastast 5_M_W_vs_8_M_W
# up
#df3[5,1] <- "M_W"
#df3[5,2] <- "Up"
#df3[5,3] <- dim(contrast_3[contrast_3$log2FoldChange > 1,])[1]
#names <- rownames(contrast_3[contrast_3$log2FoldChange > 1,])
#write.table(names, paste0(contrast_results_dir, "/","5_M_W_vs_8_M_W", "/", "5_M_W_vs_8_M_W_up"), col.names = F, row.names = F, quote = F)
# down
#df3[6,1] <- "M_W"
#df3[6,2] <- "Down"
#df3[6,3] <- dim(contrast_3[contrast_3$log2FoldChange < 1,])[1]
#names <- rownames(contrast_3[contrast_3$log2FoldChange < 1,])
#write.table(names, paste0(contrast_results_dir, "/","5_M_W_vs_8_M_W", "/", "5_M_W_vs_8_M_W_down"), col.names = F, row.names = F, quote = F)
# contrast 5_M_D_vs_5_M_W
# down
#df3[7,1] <- "5M_W-D"
#df3[7,2] <- "Down"
#df3[7,3] <- dim(contrast_4[contrast_4$log2FoldChange < 1,])[1]
#names <- rownames(contrast_4[contrast_4$log2FoldChange < 1,])
#write.table(names, paste0(contrast_results_dir, "/","5_M_D_vs_5_M_W", "/", "5_M_D_vs_5_M_W_down"), col.names = F, row.names = F, quote = F)
# up
#df3[8,1] <- "5M_W-D"
#df3[8,2] <- "Up"
#df3[8,3] <- dim(contrast_4[contrast_4$log2FoldChange > 1,])[1]
#names <- rownames(contrast_4[contrast_4$log2FoldChange > 1,])
#write.table(names, paste0(contrast_results_dir, "/","5_M_D_vs_5_M_W", "/", "5_M_D_vs_5_M_W_up"), col.names = F, row.names = F, quote = F)
# contrast 8_M_W_vs_8_M_W
# down
#df3[9,1] <- "8M_W-D"
#df3[9,2] <- "Down"
#df3[9,3] <- dim(contrast_5[contrast_5$log2FoldChange < 1,])[1]
#names <- rownames(contrast_5[contrast_5$log2FoldChange < 1,])
#write.table(names, paste0(contrast_results_dir, "/","8_M_D_vs_8_M_W", "/", "8_M_D_vs_8_M_W_down"), col.names = F, row.names = F, quote = F)
# up 
#df3[10,1] <- "8M_W-D"
#df3[10,2] <- "Up"
#df3[10,3] <- dim(contrast_5[contrast_5$log2FoldChange > 1,])[1]
#names <- rownames(contrast_5[contrast_5$log2FoldChange > 1,])
#write.table(names, paste0(contrast_results_dir, "/","8_M_D_vs_8_M_W", "/", "8_M_D_vs_8_M_W_up"), col.names = F, row.names = F, quote = F)

###
#df3$Contrast <- factor(df3$Contrasts, levels = c("10 B0", "10 B", "10 M", "10 P", "270 B0", "270 B", "270 M", "270 P"))
#levels(df3$Contrast) <- factor(df3$Contrasts, levels = c("B0", "B", "M", "T"))
#theTable$Position <- factor(theTable$Position, levels = c(...))

library(wesanderson, lib= "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(tidyverse, lib= "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
library(svglite, lib= "/Storage/data1/jorge.munoz/NRGSC.new/libraries")

#df3$Treatment <- c("D" "NR", "NR", "NR", "NR", "NR", "NR", "NR", "R", "R", "R", "R", "R", "R", "R", "R")
#cols <- wes_palette("GrandBudapest1", 2, type = "continuous")
cols <- c("#ECCBAE", "#D69C4E")
df3$Contrast <- factor(df3$Contrast, levels = c("5_M_W_vs_5_T_D", "8_M_W_vs_8_T_D"))
DEG_PLOT <- ggplot(df3, aes(x=Contrast, y=DEGs, fill=Type, label = DEGs)) +
  geom_bar(stat="identity", position = "stack", size = 3)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +	
  coord_flip() +
  xlab("Contrast") +
#  scale_x_discrete(breaks= df3$Contrast )+
  scale_fill_manual(values = cols)+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  #facet_wrap(~Treatment, drop = T, ncol = 4, nrow = 1)+
  labs(fill="")+
  theme_bw()+
  theme( text = element_text(family = "Times new roman", size=18),
     	 legend.position = "top",
	 legend.text = element_text(size = 10),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line(colour = "#030303", size = 2),
	 axis.ticks = element_line(size = 2, color="#030303"),
         axis.text.x = element_text(angle=0, size=12, color = "black"),
	 axis.text.y = element_text(angle=0, size=12, color = "black"),
         plot.title = element_text(hjust = 0.5))


ggsave(DEG_PLOT, filename = "../results/DEG_stackplot_ok.svg", device= "svg" ,units = "cm",width = 52/1.5, height = 17,dpi = 320)
