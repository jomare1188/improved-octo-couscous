## variables
wd = "/Storage/data1/jorge.munoz/DOLORES/RNAseq/code"
outdir = "/Storage/data1/jorge.munoz/DOLORES/RNAseq/code"
libraries = "/Storage/data1/jorge.munoz/NRGSC.new/libraries"
cores = 15
# set working directory
setwd(wd)
# make results dir
dir.create(outdir)
## load libraries 
library(tximport, lib.loc = libraries)
library(DESeq2, lib.loc = libraries)
library(backports, lib.loc = libraries)
library(tidyverse, lib.loc = libraries)
library(BiocParallel, lib.loc = libraries)
library(wesanderson, lib.loc = libraries)
library(tibble, lib.loc = libraries)
library(ggrepel, lib.loc = libraries)
# parallel envirorment
register(MulticoreParam(cores))
# define fuction to make DEA analysis in DESeq2
metadata="../data/metadata.csv"
tx2gen2="../data/tx2gen.csv"

## FUNCTION TO PLOT PCA USING PLOT PCA FUCTION FROM DESEQ

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
ggsave( a, filename = path ,units = "cm",width = 15*1.3, height = 15,dpi = 320)
return(a)
}


## FUNCTION TO PLOT PCA
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

#dir.create(paste0(outdir, group1, "_vs_", group2))
sample_table <-read.table(metadata, sep = ",", header = T)
tx2gene = read.table(tx2gen2, sep = ",", col.names =c("transid","geneid"))

#filter <- sample_table %>% filter( Group == group1 | Group == group2 )
# all files
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
        design = ~ Group)
        
# filter genes that have 0 expression it at least 2 columns
dim(raw)
temp <- as.data.frame(counts(raw))
logic <-(apply(temp,c(1,2), function(x){x>0}))
filter_genes <- rowSums(logic)>2
fi <- raw[filter_genes,]
dim(fi)
temp <- NULL
# normalize using vst transformation
vst <- vst(fi)
df <- assay(raw)
# Write raw counts table
write.table(df, "../results/PCA/raw_counts.csv", col.names = T, row.names = T, sep = ",", quote = F)
##################################################################
##################### MakePCA ####################################
##################################################################

vst.pca <- prcomp(x = t(assay(vst)), scale = TRUE)
var_explained <- vst.pca$sdev^2/sum(vst.pca$sdev^2)

df <- vst.pca$x %>%  as.data.frame
rownames(df) <- sample_table$sample
df <- df %>% mutate("Genotype" = as.factor(sample_table$Genotype))
df <- df %>% mutate("Group" = as.factor(sample_table$Group))
df <- df %>% mutate("Sample" = as.factor(sample_table$Sample))
df <- df %>% mutate("Treatment" = as.factor(sample_table$Treatment))
df <- df %>% mutate("State" = as.factor(sample_table$State))

# PLOT PCA 
plot_pca(file = "../results/PCA_vst_Genotype-Group.png", data = df, color = df$Genotype , shape = df$Group, title = "PCA Genotype~Group", var = var_explained, lab_color = "Genotype")
plot_pca(file = "../results/PCA_vst_Sample-Group.png", data = df, color = df$Sample , shape = df$Group, title = "PCA Sample~Group", var = var_explained, lab_color = "Sample")
plot_pca(file = "../results/PCA_vst_Treatment-Group.png", data = df, color = df$Treatment , shape = df$Group, title = "PCA Treatment~Group", var = var_explained, lab_color = "Treatment")
plot_pca(file = "../results/PCA_vst_State-Group.png", data = df, color = df$State , shape = df$Group, title = "PCA State~Group", var = var_explained, lab_color = "State")

##################################################################
##################### make aggregated PCA ########################
##################################################################

# get mean of replicates
G5M_D <- as.data.frame(assay(raw)) %>% select("21DH5_S21", "22DH5_S22", "23DH5_S23") %>% rowMeans() %>% as.data.frame()
colnames(G5M_D) <- "G5M_D"

G8M_D <- as.data.frame(assay(raw)) %>% select("21DH8_S24", "22DH8_S25", "23DH8_S26" ) %>% rowMeans() %>% as.data.frame()
colnames(G8M_D) <- "G8M_D"

G5T_D <- as.data.frame(assay(raw)) %>% select("FIDH15000_S27", "FIDH25000_S28", "FIDH35000_S29" ) %>% rowMeans() %>% as.data.frame()
colnames(G5T_D) <- "G5T_D"

G8T_D <- as.data.frame(assay(raw)) %>% select("FIDH1800_S30", "FIDH28000_S31", "FIDH38000_S32" ) %>% rowMeans() %>% as.data.frame()
colnames(G8T_D) <- "G8T_D"

G5M_W <- as.data.frame(assay(raw)) %>% select("F25C1_S6", "F25C2_S7", "F25C3_S8" ) %>% rowMeans() %>% as.data.frame()
colnames(G5M_W) <- "G5M_W"

G8M_W <- as.data.frame(assay(raw)) %>% select("F28C1_S9", "F28C2_S10", "F28C3_S11" ) %>% rowMeans() %>% as.data.frame()
colnames(G8M_W) <- "G8M_W"
# make new matrix
aggregated <- cbind(G5M_D, G8M_D, G5T_D, G8T_D, G5M_W, G8M_W)

# write aggregated raw matrix
write.table(aggregated, "../results/PCA/raw_aggragated_counts.csv", col.names = T, row.names = T, sep = ",", quote = F)

#temp <- as.data.frame(counts(aggregated))
logic <-(apply(aggregated,c(1,2), function(x){x>0}))
filter_genes <- rowSums(logic)>1
fi <- round(as.matrix(aggregated[filter_genes,]))
dim(fi)
# normalize using vst transformation
vst <- vst(fi)

# PCA
vst.pca <- prcomp(x = t(vst))
var_explained <- vst.pca$sdev^2/sum(vst.pca$sdev^2)
df <- vst.pca$x %>%  as.data.frame

#rownames(df) <- sample_table$sample
df <- df %>% mutate("Genotype" = as.factor(c("5000", "8008", "5000", "8008", "5000", "8008")))
df <- df %>% mutate("Group" = as.factor(c("G5M_D", "G8M_D", "G5T_D", "G8T_D", "G5M_W", "G8M_W")))
df <- df %>% mutate("Treatment" = as.factor(c("Drought", "Drought", "Drought", "Drought", "Watered", "Watered")))
df <- df %>% mutate("State" = as.factor(c("Maturation", "Maturation", "Tillering", "Tillering", "Maturation", "Maturation")))

#PLOT PCA
plot_pca(file = "../results/Agreggated_PCA_vst_Genotype-Group.png", data = df, color = df$Genotype , shape = df$Group, title = "PCA Genotype~Group", var = var_explained, lab_color = "Genotype")
plot_pca(file = "../results/Agreggated_PCA_vst_Treatment-Group.png", data = df, color = df$Treatment , shape = df$Group, title = "PCA Treatment~Group", var = var_explained, lab_color = "Treatment")
plot_pca(file = "../results/Agreggated_PCA_vst_State-Group.png", data = df, color = df$State , shape = df$Group, title = "PCA State~Group", var = var_explained, lab_color = "State")

# PLOT PCA USING DESEQ FUNCTION
plot_pca_deseq(data = vst, shape = vst$Group, title = "PCA Genotype~Group", name = "Genotype", path = "../results/DESEQ_PCA_vst_Genotype-Group.png")
plot_pca_deseq(data = vst, shape = vst$Group, title = "PCA State~Group", name = "State", path = "../results/DESEQ_PCA_vst_State-Group.png")
plot_pca_deseq(data = vst, shape = vst$Group, title = "PCA Treatment~Group", name = "Treatment", path = "../results/DESEQ_PCA_vst_Treatment-Group.png")
plot_pca_deseq(data = vst, shape = vst$Group, title = "PCA Sample~Group", name = "Sample", path = "../results/DESEQ_PCA_vst_Sample-Group.png")
