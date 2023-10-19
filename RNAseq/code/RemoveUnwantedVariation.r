wd = "/Storage/data1/jorge.munoz/DOLORES/RNAseq/code"
outdir = "/Storage/data1/jorge.munoz/DOLORES/RNAseq/results"
libraries = "/Storage/data1/jorge.munoz/NRGSC.new/libraries"

# set working directory
setwd(wd)
# make results dir
dir.create(outdir)
## load libraries 

library(readr, lib.loc =  libraries)
library(dplyr, lib.loc =  libraries)
library(magrittr,  lib.loc = libraries)
library(tximport, lib.loc = libraries)
library(DESeq2, lib.loc = libraries)
library(ggplot2, lib.loc =  libraries)
library(hexbin, lib.loc =  libraries)
library(EDASeq, lib.loc =  libraries)
library(RUVSeq, lib.loc =  libraries)

metadata="../data/metadata.csv"
tx2gen2="../data/tx2gen.csv"

plot_pca_deseq <- function(data, shape, title, name, path) {
a <- plotPCA(data, intgroup = name) +
theme_bw() +
geom_point(size=4.5, aes( shape = shape))+
scale_shape_manual(values = seq(0,8)) +
labs(title = title, col= name, shape = "Group") +
theme( text = element_text(size=22),
panel.border = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
plot.title = element_text(hjust = 0.5),
axis.line = element_line(colour = "black") )
ggsave( a, filename = path ,units = "cm",width = 15*1.3, height = 15*1.4,dpi = 320)
return(a)
}


sample_table <-read.table(metadata, sep = ",", header = T)
tx2gene = read.table(tx2gen2, sep = ",", col.names =c("transid","geneid"))

sample_files = paste0("/Storage/data1/jorge.munoz/DOLORES/quant_test/longest_transcript_orthogroup/results/", pull(sample_table, "Files"), "/quant.sf")
# filter by groups
# sample_files = paste0("/Storage/data1/jorge.munoz/NRGSC/", pull(filter , "Sample_file"), "/quant.sf")
# name table columns
names(sample_files) = pull(sample_table, "Files")
count_data = tximport( files = sample_files,
        type = "salmon",
        tx2gene =  tx2gene,
        ignoreTxVersion = F)

sample_table$Treatment <- as.factor((sample_table$Treatment))
sample_table$Genotype <- as.factor((sample_table$Genotype))
sample_table$State <- as.factor((sample_table$State))
sample_table$Replicate <- as.factor((sample_table$Replicate))
sample_table$Sample <- as.factor((sample_table$Sample))
sample_table$Group <- as.factor((sample_table$Group))

raw <- DESeqDataSetFromTximport(txi = count_data,
        colData = sample_table,
        design = ~ Genotype)

dim(raw)

# RUVg 
EC = 1000
dea <- DESeq(raw)

results_table <- as.data.frame(results(dea))
filtered_results_table <- results_table %>% filter(padj < 0.05 ) %>% arrange(abs(log2FoldChange))

empirical_control <- rownames(head(filtered_results_table), EC)
RUVgSet <- RUVg(as.matrix(assay(raw)), empirical_control, k = 1)
dataRUVg<- DESeqDataSetFromMatrix(countData = RUVgSet$normalizedCounts,
                                 colData = sample_table,
                                 design = ~ Genotype)

dataRUVg <- estimateSizeFactors(dataRUVg)
RUVg_vst <- vst(dataRUVg)

plot_pca_deseq(data = RUVg_vst, shape = RUVg_vst$Group, title = "RUVg PCA Genotype~Group", name = "Genotype", path = "../results/RUVg_Genotype-Group.png")

## RUVs

differences <- makeGroups(sample_table$Group)
genes <- row.names(raw)
setRUVs <- RUVs(as.matrix(assay(raw)), genes, k=1, differences)
dataRUVs <- DESeqDataSetFromMatrix(countData = setRUVs$normalizedCounts,
                                  colData = sample_table,
                                  design = ~ Group)
dataRUVs <- estimateSizeFactors(dataRUVs)
RUVs_vst <- vst(dataRUVs)

plot_pca_deseq(data = RUVs_vst, shape = RUVs_vst$Group, title = "RUVs PCA Genotype~Group", name = "Genotype", path = "../results/RUVs_Genotype-Group.png")

# RUVr residuals

design <- model.matrix(~sample_table$Group)
y <- DGEList(counts=counts(raw), group=sample_table$Group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

set_RUVr <- RUVr(y$counts,rownames(y), k=1, res)
dataRUVr <- DESeqDataSetFromMatrix(countData = set_RUVr$normalizedCounts,
                                   colData = sample_table,
                                   design = ~ Group)

RUVr_vst <- varianceStabilizingTransformation(dataRUVr)
plot_pca_deseq(data = RUVr_vst, shape = RUVr_vst$Group, title = "RUVr PCA Genotype~Group", name = "Genotype", path = "../results/RUVr_Genotype-Group.png")

