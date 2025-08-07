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

# Extract and save TPM counts
tpm_counts <- count_data$abundance
write.table(tpm_counts, paste0(outdir, "tpm_counts.csv"), 
            col.names = TRUE, row.names = TRUE, sep = ",", quote = FALSE)

# Prepare sample metadata
sample_table$Treatment <- as.factor(sample_table$Treatment)
sample_table$Genotype <- as.factor(sample_table$Genotype)
sample_table$Group <- as.factor(sample_table$Group)
sample_table$Generation <- as.factor(sample_table$Generation)
sample_table$Tissue <- as.factor(sample_table$Tissue)

# DESeq2 Analysis ---------------------------------------------------------

# Create DESeq2 object
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

# Variance stabilizing transformation
vst_all <- vst(fi)

# Save count matrices
df <- assay(raw)
write.table(df, paste0(outdir, "raw_counts.csv"), 
            col.names = TRUE, row.names = TRUE, sep = ",", quote = FALSE)

tpm_filtered <- tpm_counts[filter_genes, ]
write.table(tpm_filtered, paste0(outdir, "tpm_counts_filtered.csv"), 
            col.names = TRUE, row.names = TRUE, sep = ",", quote = FALSE)
write.table(tpm_counts, paste0(outdir, "tpm_counts_raw.csv"), 
            col.names = TRUE, row.names = TRUE, sep = ",", quote = FALSE)

# Create tissue-specific VST matrices
vst_mothers <- vst(fi[, fi$Tissue == "Bud"])
vst_sons <- vst(fi[, fi$Tissue == "Leaf"])

# Save tissue-specific VST matrices
write.table(assay(vst_mothers), paste0(outdir, "vst_mothers.csv"), 
            col.names = TRUE, row.names = TRUE, sep = ",", quote = FALSE)
write.table(assay(vst_sons), paste0(outdir, "vst_sons.csv"), 
            col.names = TRUE, row.names = TRUE, sep = ",", quote = FALSE)

# PCA Analysis ------------------------------------------------------------

cat("Generating PCA plots...\n")

# PCA plots with all samples
plot_pca_deseq_custom(data = vst_all, 
                      title = "PCA Genotype - Group - Treatment",
                      grouping_vars = c("Group", "Genotype", "Treatment"),
                      n_top = 500,
                      color_var = "Group",
                      shape_var = "Treatment",
                      size_var = "Genotype",
                      path = paste0(outdir, "DESEQ_PCA_vst_Genotype-Group-Treatment.png"))

plot_pca_deseq_custom(data = vst_all, 
                      title = "PCA Generation - Treatment - Genotype",
                      grouping_vars = c("Generation", "Treatment", "Genotype"),
                      n_top = 500,
                      color_var = "Treatment",
                      shape_var = "Generation",
                      size_var = "Genotype",
                      path = paste0(outdir, "DESEQ_PCA_vst_Treatment-Generation-Genotype.png"))

# PCA plots for mothers (Bud tissue)
plot_pca_deseq_custom(data = vst_mothers, 
                      title = "PCA Mothers: Genotype - Group - Treatment",
                      grouping_vars = c("Group", "Genotype", "Treatment"),
                      n_top = 500,
                      color_var = "Group",
                      shape_var = "Treatment",
                      size_var = "Genotype",
                      path = paste0(outdir, "mothers_DESEQ_PCA_vst_Genotype-Group-Treatment.png"))

plot_pca_deseq_custom(data = vst_mothers, 
                      title = "PCA Mothers: Phenological stage - Treatment - Genotype",
                      grouping_vars = c("Phenological_stage", "Treatment", "Genotype"),
                      n_top = 500,
                      color_var = "Treatment",
                      shape_var = "Phenological_stage",
                      size_var = "Genotype",
                      path = paste0(outdir, "mothers_DESEQ_PCA_vst_Treatment-Phenological_stage-Genotype.png"))

# PCA plots for sons (Leaf tissue)
plot_pca_deseq_custom(data = vst_sons, 
                      title = "PCA Sons: Genotype - Group - Treatment",
                      grouping_vars = c("Group", "Genotype", "Treatment"),
                      n_top = 500,
                      color_var = "Group",
                      shape_var = "Treatment",
                      size_var = "Genotype",
                      path = paste0(outdir, "sons_DESEQ_PCA_vst_Genotype-Group-Treatment.png"))

plot_pca_deseq_custom(data = vst_sons, 
                      title = "PCA Sons: Generation - Treatment - Genotype",
                      grouping_vars = c("Generation", "Treatment", "Genotype"),
                      n_top = 500,
                      color_var = "Generation",
                      shape_var = "Treatment",
                      size_var = "Genotype",
                      path = paste0(outdir, "sons_DESEQ_PCA_vst_Generation-Treatment-Genotype.png"))

cat("Analysis complete! Check", outdir, "for results.\n")
