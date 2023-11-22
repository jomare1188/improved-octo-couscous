## variables
wd = "/Storage/data1/jorge.munoz/DOLORES/RNAseq/code"
outdir = "../results/"
cores = 6
# set working directory
setwd(wd)
# make results dir
dir.create(outdir)
libraries = "/Storage/data1/jorge.munoz/NRGSC.new/libraries"
## load libraries 
library(tximport, lib.loc = libraries)
library(DESeq2, lib.loc = libraries)
library(backports, lib.loc = libraries)
library(dplyr, lib.loc = libraries)
library(tidyverse, lib.loc = libraries)
library(BiocParallel, lib.loc = libraries)

# parallel envirorment
register(MulticoreParam(cores))

# define fuction to make DEA analysis in DESeq2
dea_by_group <- function(group1, group2) {
	dir.create(paste0(outdir, group1, "_vs_", group2))
        sample_table <-read.table("../data/metadata.csv", sep = ",", header = T)      
        filter <- sample_table %>% filter( Group == group1 | Group == group2 )
	
        # load files paths
	sample_files = paste0("/Storage/data1/jorge.munoz/DOLORES/quant_test/longest_transcript_orthogroup/results/", pull(filter, "Files"), "/quant.sf")
        # name table columns
        names(sample_files) = pull(filter, "Files")
        # relate genes to transcripts
        tx2gene = read.table("../data/tx2gen.csv", sep = ",", col.names =c("transid","geneid"))
	# GENE MODE 
        # import count data to tximport
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
                                        colData = filter,
                                        design = ~ Group)
        dim(raw)
        temp <- as.data.frame(counts(raw))
        logic <-(apply(temp,c(1,2), function(x){x>0}))
        filter_genes <- rowSums(logic)>2
        fi <- raw[filter_genes,]
        dim(fi)
        temp <- NULL
        data <- estimateSizeFactors(fi)
        vst <- varianceStabilizingTransformation(data)
        
        ### Differencial expression analyses
        dea <- DESeq(fi, parallel = T)
        dea_contrast <- results(dea, lfcThreshold= 1, altHypothesis="greaterAbs", parallel = T)
        #save(dea_contrast, file = paste0("./DEA/by_groups/", group1,"_vs_", group2, ".Rdata", sep = ""))
        
	dea_df <- as.data.frame(dea_contrast)
	##  Save files for other things

	baseMeanA <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$Group == group1])
        baseMeanB <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$Group == group2])
        
        res = cbind(baseMeanA, baseMeanB, dea_df)
        res = cbind(sampleA=group1, sampleB=group2, as.data.frame(res))
	res = res[complete.cases(res),]
        res = as.data.frame(res[order(res$pvalue),])
	
	filtered <- res %>% filter(pvalue < 0.05 & abs(log2FoldChange) > 1)

        write.table(filtered, file = paste0(outdir, "/",  group1, "_vs_", group2, "/" ,group1, "_vs_", group2, ".DESeq.DE_results"), sep = "\t", quote = F )
        ### More tables
        df_fi <- assay(fi)
        df_vst <- assay(vst)	
	write.table(df_fi, file = paste0(outdir, "/",  group1, "_vs_", group2, "/", group1, "-", group2, "-", "filtered_raw_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
	write.table(df_vst, file = paste0(outdir, "/",  group1, "_vs_", group2,"/" , group1, "-", group2, "-", "vst_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
	write.table(filter, file = paste0(outdir, "/",  group1, "_vs_", group2,"/", group1, "-", group2, "-", "metadata.tsv", sep = ""), sep = "\t", quote = F, col.names = T)
}

dea_by_group("5_M_D", "8_M_D")
dea_by_group("5_M_W", "8_M_W")
dea_by_group("5_T_D", "8_T_D")
dea_by_group("5_M_D", "5_M_W")
dea_by_group("8_M_D", "8_M_W")
# NEW CONTRASTAS
dea_by_group("5_M_W", "5_T_D")
dea_by_group("8_M_W", "8_T_D")




