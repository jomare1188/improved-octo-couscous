## variables
cores = 5
# make results dir
dir.create(outdir)
## load libraries 
library(tximport)
library(DESeq2)
library(tidyverse)
library(BiocParallel)
library(wesanderson)
library(tibble)
library(ggrepel)
# parallel envirorment
register(MulticoreParam(cores))
metadata="../data/ALL_meta_data_dolores.csv"
tx2gen2="../data/tx2gen_longest.csv"


## Function to plot PCA from raw R
plot_pca <- function(data, color, shape, title, var, file, lab_color){
ggplot(df, aes(x=PC1, y=PC2, colour = color, shape = shape)) +
geom_point(size = 4.5) +
theme_bw(base_size=20) +
labs(title = title, shape = "Group", color = lab_color) +
labs(x=paste0("PC1: ",round(var[1]*100,1),"%"),
        y=paste0("PC2: ",round(var[2]*100,1),"%")) +
theme( text = element_text(size=20),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
        axis.text.y=element_text(colour="black"),
        axis.line = element_line(colour = "black"))

ggsave(filename = file ,units = "cm", width = 25, height = 25,dpi = 320)
}

# Alternative version that lets you specify which variables map to which aesthetics using deseq2 function plotPCA
plot_pca_deseq_custom <- function(data, title, grouping_vars, n_top, path,
                                  color_var = NULL, shape_var = NULL,
                                  size_var = NULL) {

  pca_data <- plotPCA(data, intgroup = grouping_vars, returnData = TRUE, ntop = n_top)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  # Create aesthetic mapping
  aes_mapping <- aes(PC1, PC2)

  if (!is.null(color_var) && color_var %in% grouping_vars) {
    aes_mapping$colour <- sym(color_var)
    aes_mapping$fill <- sym(color_var)
  }

  if (!is.null(shape_var) && shape_var %in% grouping_vars) {
    aes_mapping$shape <- sym(shape_var)
  }

  if (!is.null(size_var) && size_var %in% grouping_vars) {
    aes_mapping$size <- sym(size_var)
  }

  # Base plot with dynamic aesthetics
  p <- ggplot(pca_data, aes_mapping) +
    theme_bw() +
    labs(title = title,
         x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 22),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black"))

  # Add geom_point with conditional size
  if (!is.null(size_var) && size_var %in% grouping_vars) {
    p <- p + geom_point(stroke = 1)  # No fixed size when using size aesthetic
  } else {
    p <- p + geom_point(size = 4.5, stroke = 1)  # Fixed size when not using size aesthetic
  }

  # Add labels and scales
  if (!is.null(color_var) && color_var %in% grouping_vars) {
    p <- p + labs(color = color_var) +
      guides(fill = "none")
  }

  if (!is.null(shape_var) && shape_var %in% grouping_vars) {
    p <- p + labs(shape = shape_var) +
      scale_shape_manual(values = c(21, 22, 23, 24, 25, 8))
  }

  if (!is.null(size_var) && size_var %in% grouping_vars) {
    p <- p + labs(size = size_var) +
      scale_size_manual(values = c(4, 6))  # Customize sizes for better differentiation
  }

  # Save plot
  ggsave(plot = p, filename = path, units = "cm",
         width = 15*2, height = 15*2, dpi = 320)

  return(p)
}

#dir.create(paste0(outdir, group1, "_vs_", group2))
sample_table <-read.table(metadata, sep = ",", header = T) %>% filter(Exist == "YES")
tx2gene = read.table(tx2gen2, sep = ",", col.names =c("transid","geneid"))

#filter <- sample_table %>% filter( Group == group1 | Group == group2 )
# all files
sample_files = paste0("/Storage/data1/jorge.munoz/DOLORES/nf/R2C/results/9_salmonQuant/", pull(sample_table, "Sample"), "/quant.sf")
# filter by groups
# sample_files = paste0("/Storage/data1/jorge.munoz/NRGSC/", pull(filter , "Sample_file"), "/quant.sf")
# name table columns

names(sample_files) = pull(sample_table, "Sample")
count_data = tximport( files = sample_files, type = "salmon", tx2gene =  tx2gene, ignoreTxVersion = F)
        
# Extract TPM counts from tximport object
tpm_counts <- count_data$abundance


# Write TPM counts table
write.table(tpm_counts, "../results/PCA/tpm_counts.csv", col.names = T, row.names = T, sep = ",", quote = F)


sample_table$Treatment <- as.factor((sample_table$Treatment))
sample_table$Genotype <- as.factor((sample_table$Genotype))
#sample_table$Phenological_stage <- as.factor((sample_table$State))
#sample_table$Repetition <- as.factor((sample_table$Replicate))
#sample_table$Sample <- as.factor((sample_table$Sample))
sample_table$Group <- as.factor((sample_table$Group))
sample_table$Generation <- as.factor((sample_table$Generation))
sample_table$Tissue <- as.factor((sample_table$Tissue))



raw <- DESeqDataSetFromTximport(txi = count_data,
	colData = sample_table,
        design = ~ Group)
        
# filter genes that have 0 expression it at least 2 columns
dim(raw)
temp <- as.data.frame(counts(raw))
logic <-(apply(temp,c(1,2), function(x){x>0}))
filter_genes <- rowSums(logic)>5
fi <- raw[filter_genes,]
dim(fi)
temp <- NULL

# normalize using vst transformation
vst <- vst(fi)
df <- assay(raw)

# Write raw counts table
write.table(df, "../results/PCA/raw_counts.csv", col.names = T, row.names = T, sep = ",", quote = F)

# Write tmp counts
tpm_filtered <- tpm_counts[filter_genes,]
write.table(tpm_filtered, "../results/PCA/tpm_counts_filtered.csv", col.names = T, row.names = T, sep = ",", quote = F)
write.table(tpm_counts, "../results/PCA/tpm_counts_raw.csv", col.names = T, row.names = T, sep = ",", quote = F)

# filter vst matrices for each batch 
vst_mothers <- vst(fi[, fi$Tissue == "Bud"])
vst_sons <- vst(fi[, fi$Tissue == "Leaf"])

# write vst matrices for each batch
write.table(assay(vst_mothers), "../results/PCA/vst_mothers.csv", col.names = T, row.names = T, sep = ",", quote = F)
write.table(assay(vst_sons), "../results/PCA/vst_sons.csv", col.names = T, row.names = T, sep = ",", quote = F)



# plot PCAs with all samples together
plot_pca_deseq_custom(data = vst, 
                      title = "PCA Genotype - Group - Treatment",
                      grouping_vars = c("Group", "Genotype", "Treatment"),
		      n_top = 500,
                      color_var = "Group",
                      shape_var = "Treatment",
		      size_var = "Genotype",
                      path = "../results/PCA/DESEQ_PCA_vst_Genotype-Group-Treatment.png")

plot_pca_deseq_custom(data = vst, 
                      title = "PCA Generation - Treatment - Genotype",
                      grouping_vars = c("Generation", "Treatment", "Genotype"),
		      n_top = 500,
                      color_var = "Treatment",
                      shape_var = "Generation",
		      size_var = "Genotype",
                      path = "../results/PCA/DESEQ_PCA_vst_Treatment-Generation-Genotype.png")

# plot PCAs for each batch of experiments (mothers and sons)
# mothers
plot_pca_deseq_custom(data = vst_mothers, 
                      title = "PCA Genotype - Group - Treatment",
                      grouping_vars = c("Group", "Genotype", "Treatment"),
		      n_top = 500,
                      color_var = "Group",
                      shape_var = "Treatment",
		      size_var = "Genotype",
                      path = "../results/PCA/mothers_DESEQ_PCA_vst_Genotype-Group-Treatment.png")

plot_pca_deseq_custom(data = vst_mothers, 
                      title = "PCA Phenological stage - Treatment - Genotype",
                      grouping_vars = c("Phenological_stage", "Treatment", "Genotype"),
		      n_top = 500,
                      color_var = "Treatment",
                      shape_var = "Phenological_stage",
		      size_var = "Genotype",
                      path = "../results/PCA/mothers_DESEQ_PCA_vst_Treatment-Phenological_stage-Genotype.png")

# sons 
plot_pca_deseq_custom(data = vst_sons, 
                      title = "PCA Genotype - Group - Treatment",
                      grouping_vars = c("Group", "Genotype", "Treatment"),
		      n_top = 500,
                      color_var = "Group",
                      shape_var = "Treatment",
		      size_var = "Genotype",
                      path = "../results/PCA/sons_DESEQ_PCA_vst_Genotype-Group-Treatment.png")

plot_pca_deseq_custom(data = vst_sons, 
                      title = "PCA Generation - Treatment - Genotype",
                      grouping_vars = c("Generation", "Treatment", "Genotype"),
		      n_top = 500,
                      color_var = "Generation",
                      shape_var = "Treatment",
		      size_var = "Genotype",
                      path = "../results/PCA/sons_DESEQ_PCA_vst_Generation-Treatment-Genotype.png")


