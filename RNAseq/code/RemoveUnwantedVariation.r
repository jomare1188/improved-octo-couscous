# RNA-seq Preliminary Analysis Script
# DOLORES Project
# Author: Your Name

# Setup -------------------------------------------------------------------

## Variables
cores <- 5
outdir <- "../results/PCA/"  # Define your output directory

# Create results directory
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Load libraries
library(tximport)
library(DESeq2)
library(tidyverse)
library(BiocParallel)
library(wesanderson)
library(tibble)
library(ggrepel)
library(RUVSeq)

# Set up parallel environment
register(MulticoreParam(cores))

# Source custom functions
source("pca_functions.r")

# File paths
metadata <- "../data/ALL_meta_data_dolores.csv"
tx2gen2 <- "../data/tx2gen_longest.csv"

# Data Import and Processing ----------------------------------------------

# Load metadata and gene mapping
sample_table <- read.table(metadata, sep = ",", header = TRUE) %>%
  filter(Exist == "YES")

tx2gene <- read.table(tx2gen2, sep = ",", col.names = c("transid", "geneid"))

# Prepare file paths for salmon quantification files
sample_files <- paste0("/Storage/data1/jorge.munoz/DOLORES/nf/R2C/results/9_salmonQuant/",
                       pull(sample_table, "Sample"), "/quant.sf")
names(sample_files) <- pull(sample_table, "Sample")

# Import count data using tximport
count_data <- tximport(files = sample_files,
                       type = "salmon",
                       tx2gene = tx2gene,
                       ignoreTxVersion = FALSE)

# Prepare sample metadata
sample_table$Treatment <- as.factor(sample_table$Treatment)
sample_table$Genotype <- as.factor(sample_table$Genotype)
sample_table$Group <- as.factor(sample_table$Group)
sample_table$Generation <- as.factor(sample_table$Generation)
sample_table$Tissue <- as.factor(sample_table$Tissue)

# DESeq2 Analysis ---------------------------------------------------------

raw <- DESeqDataSetFromTximport(txi = count_data,
                                colData = sample_table,
                                design = ~ Group)

# Quality control: filter genes with low expression
cat("Original dimensions:", dim(raw), "\n")
temp <- as.data.frame(counts(raw))
logic <- apply(temp, c(1, 2), function(x) {x > 0})
filter_genes <- rowSums(logic) > 5
fi <- raw[filter_genes, ]
cat("Filtered dimensions:", dim(fi), "\n")
rm(temp)  # Clean up memory


# RUVg ---------------------------------------------------------------------
EC = 1000
dea <- DESeq(fi)

results_table <- as.data.frame(results(dea))
filtered_results_table <- results_table %>% filter(padj < 0.05 ) %>% arrange(abs(log2FoldChange))

empirical_control <- rownames(head(filtered_results_table), EC)
RUVgSet <- RUVg(as.matrix(assay(fi)), empirical_control, k = 1)
dataRUVg<- DESeqDataSetFromMatrix(countData = RUVgSet$normalizedCounts,
                                 colData = sample_table,
                                 design = ~ Genotype)

dataRUVg <- estimateSizeFactors(dataRUVg)
RUVg_vst <- vst(dataRUVg)

# plot PCAs with all samples together
plot_pca_deseq_custom(data = RUVg_vst, 
                      title = "PCA Genotype - Group - Treatment",
                      grouping_vars = c("Group", "Genotype", "Treatment"),
		      n_top = 500,
                      color_var = "Group",
                      shape_var = "Treatment",
		      size_var = "Genotype",
                      path = "../results/PCA/RUVg_Genotype-Group-Treatment.png")

plot_pca_deseq_custom(data = RUVg_vst, 
                      title = "PCA Generation - Treatment - Genotype",
                      grouping_vars = c("Generation", "Treatment", "Genotype"),
		      n_top = 500,
                      color_var = "Treatment",
                      shape_var = "Generation",
		      size_var = "Genotype",
                      path = "../results/PCA/RUVg_Treatment-Generation-Genotype.png")

## RUVs --------------------------------------------------------------------------

differences <- makeGroups(sample_table$Group)
genes <- row.names(fi)
setRUVs <- RUVs(as.matrix(assay(fi)), genes, k=1, differences)
dataRUVs <- DESeqDataSetFromMatrix(countData = setRUVs$normalizedCounts,
                                  colData = sample_table,
                                  design = ~ Group)
dataRUVs <- estimateSizeFactors(dataRUVs)
RUVs_vst <- vst(dataRUVs)

# plot PCAs with all samples together
plot_pca_deseq_custom(data = RUVs_vst, 
                      title = "PCA Genotype - Group - Treatment",
                      grouping_vars = c("Group", "Genotype", "Treatment"),
		      n_top = 500,
                      color_var = "Group",
                      shape_var = "Treatment",
		      size_var = "Genotype",
                      path = "../results/PCA/RUVs_Genotype-Group-Treatment.png")

plot_pca_deseq_custom(data = RUVs_vst, 
                      title = "PCA Generation - Treatment - Genotype",
                      grouping_vars = c("Generation", "Treatment", "Genotype"),
		      n_top = 500,
                      color_var = "Treatment",
                      shape_var = "Generation",
		      size_var = "Genotype",
                      path = "../results/PCA/RUVs_Treatment-Generation-Genotype.png")

# RUVr residuals --------------------------------------------------------------------

design <- model.matrix(~sample_table$Group)
y <- DGEList(counts=counts(fi), group=sample_table$Group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

# Loop several k values

for (i in 1:6) {
  # Create RUVr set with current k value
  set_RUVr <- RUVr(y$counts, rownames(y), k=i, res)

  # Create DESeq dataset
  dataRUVr <- DESeqDataSetFromMatrix(countData = set_RUVr$normalizedCounts,
                                     colData = sample_table,
                                     design = ~ Group)

  # Apply variance stabilizing transformation
  RUVr_vst <- varianceStabilizingTransformation(dataRUVr)

  # Generate PCA plot with k value in filename
  plot_pca_deseq_custom(data = RUVr_vst,
	         title = paste0("PCA Generation - Treatment - Genotype", " k = ", i),
		 grouping_vars = c("Generation", "Treatment", "Genotype"),
		 n_top = 500,
		 color_var = "Treatment",
		 shape_var = "Generation",
                 size_var = "Genotype",
                 path = paste0("../results/PCA/RUVr_Treatment-Generation-Genotype_k", i, ".png"))

  plot_pca_deseq_custom(data = RUVr_vst,
                 title = paste0("PCA Genotype - Group - Treatment", " k = ", i),
		 grouping_vars = c("Genotype", "Group",  "Treatment"),
		 n_top = 500,
		 color_var = "Group",
		 shape_var = "Treatment",
                 size_var = "Genotype",
                 path = paste0("../results/PCA/RUVr_Genotype-Group-Treatment_k", i, ".png"))
}
