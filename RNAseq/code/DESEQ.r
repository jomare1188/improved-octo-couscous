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

## FUNCTION TO PLOT PCA
plot_pca <- function(data, color, shape, title, var, file, lab_color){
ggplot(df, aes(x=PC1, y=PC2, colour = color, shape = shape)) +
geom_point(size = 4.5) +
theme_bw(base_size=20) +
labs(title = "PCA Genotype~Group", shape = "Group", color = lab_color) +
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

dir.create(paste0(outdir, group1, "_vs_", group2))
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

##################################################################
##################### Differential expression analysis ###########
##################################################################

  

        ### Differencial expression analyses
        dea <- DESeq(fi, parallel = T)
        dea_contrast <- results(dea, lfcThreshold= 1, altHypothesis="greaterAbs", parallel = T)
        #save(dea_contrast, file = paste0("./DEA/by_groups/", group1,"_vs_", group2, ".Rdata", sep = ""))
        
	dea_df <- as.data.frame(dea_contrast)
	##  Save files for trinotate
	
	baseMeanA <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$Group == group1])
        baseMeanB <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$Group == group2])
        
        res = cbind(baseMeanA, baseMeanB, dea_df)
        res = cbind(sampleA=group1, sampleB=group2, as.data.frame(res))
	res = res[complete.cases(res),]
        res = as.data.frame(res[order(res$pvalue),])
        write.table(res, file = paste0(outdir, group1, "_vs_", group2,"/Trinity.isoform.counts.matrix.", group1, "_vs_", group2, ".DESeq.DE_results"), sep = "\t", quote = F )
	write.table(res, file = paste0(outdir,"Trinity.isoform.counts.matrix.", group1, "_vs_", group2, ".DESeq.DE_results"), sep = "\t", quote = F )

        write.table(cbind(as.character(filter$Group), filter$Sample_file), file= paste0(outdir, group1, "_vs_", group2,"/Trinity.isoform.counts.matrix.", group1, "_vs_", group2, ".DESeq.samples"), col.names = F, row.names = F, sep = "\t", quote = F)
        
        ### filter DEGs 
        filter_1 <- dea_df[complete.cases(dea_df),]
        filter_2 <- filter_1[filter_1$padj < 0.05,]
        filter_3 <- filter_2[abs(filter_2$log2FoldChange) > 1,]
	
        ## PlotMA  
        ylim <- c(-10, 10)
        png(paste0(outdir,"plotMA-", group1, "-", group2, ".png", sep=""),  width = 15*1.3, height = 15, res = 320, units = "cm", pointsize = 12, bg = "white")
        DESeq2::plotMA(dea_contrast, ylim=ylim, main = paste0( group1 ,"-", group2, sep = "" ), alpha = 0.05) 
        dev.off()
        ## volcano plot
        filter_1$test = filter_1$padj < 0.05 & abs(filter_1$log2FoldChange) > 1
        filter_1 = rownames_to_column(filter_1, var='ensgene')
        
        g = ggplot(filter_1, aes(x = log2FoldChange,
                                 y = -log10(padj), 
                                 name = ensgene)) +
                geom_point(aes(colour = test), size = 1, alpha = 0.3) +
                scale_colour_manual(values=c('black', 'red')) +
                geom_vline(xintercept = 1, colour ='green', linetype = 3) +
                geom_vline(xintercept = -1, colour ='green', linetype = 3) +
                geom_hline(yintercept = -log10(0.05), colour = 'blue', linetype = 3) +
                labs(title = paste( group1, "vs", group2)) +
                theme_bw() +
                theme( text = element_text(size = 22), 
                       panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"),
                       legend.position = 'none')
        
        ggsave( g, filename = paste0(outdir,"volcanoplot-",  group1, "-", group2, ".png"), units = "cm", width = 15*1.3, height = 15, dpi = 320)
        
	filter_3$names <- row.names(filter_3)
	write.table(filter_3, file = paste0(outdir, group1, "_vs_", group2, "/", group1, "-", group2, "-", "DEA.tsv", sep = ""), sep = "\t", quote = F, col.names = T)
	write.table(filter_3, file = paste0(outdir, group1, "-", group2, "-", "DEA.tsv", sep = ""), sep = "\t", quote = F, col.names = T)

        df_fi <- assay(fi)
        df_vst <- assay(vst)	

	write.table(df_fi, file = paste0(outdir, group1, "_vs_", group2, "/", group1, "-", group2, "-", "filtered_raw_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
	write.table(df_fi, file = paste0(outdir, group1, "-", group2, "-", "filtered_raw_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)

        write.table(df_vst, file = paste0(outdir, group1, "_vs_", group2, "/", group1, "-", group2, "-", "vst_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
	write.table(df_vst, file = paste0(outdir, group1, "-", group2, "-", "vst_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)

	write.table(filter, file = paste0(outdir, group1, "_vs_", group2, "/", group1, "-", group2, "-", "metadata.tsv", sep = ""), sep = "\t", quote = F, col.names = T)
	write.table(filter, file = paste0(outdir, group1, "-", group2, "-", "metadata.tsv", sep = ""), sep = "\t", quote = F, col.names = T)
}

#dea_by_group("./../data/metadata_complete.csv", "R_270_B0", "R_10_B0", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "NR_270_B0", "NR_10_B0", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "R_270_B", "R_10_B", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "NR_270_B", "NR_10_B", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "R_270_M", "R_10_M", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "NR_270_M", "NR_10_M", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "R_270_P", "R_10_P", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#dea_by_group("./../data/metadata_complete.csv", "NR_270_P", "NR_10_P", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")

dea_by_group("./../data/metadata_complete.csv", "R_270_B0", "NR_270_B0", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_10_B0", "NR_10_B0", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_270_B", "NR_270_B", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_10_B", "NR_10_B", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_270_M", "NR_270_M", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_10_M", "NR_10_M", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_270_P", "NR_270_P", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
dea_by_group("./../data/metadata_complete.csv", "R_10_P", "NR_10_P", "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")



