setwd("/Storage/data1/jorge.munoz/DOLORES/quant_test/pretty_results")
reads <- read.table("all_results.csv", sep = ",", header=FALSE)
colnames(reads) <- c('Sample', 'Category', 'Mapping_rate')
head(reads)

samples<-sort(unique(reads$Sample))


library(ggplot2, lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
counter <- 1
for (i in 1:18){
  if(i%%4==0){
    start <- i-3
    end <- i
    sams <- samples[start:end]
    ggplot(reads[which(reads$Sample %in% sams),],aes(y=Mapping_rate,x=Category,fill=Category))+
      geom_bar(stat = 'identity')+
      scale_y_sqrt()+
      facet_grid(~Sample)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=60, hjust=1))
    outfile=paste('quant_test',counter,'.pdf',sep='')
    ggsave(filename = outfile, plot=last_plot(),device = 'pdf',width = 14,height = 12)
    counter=counter+1
  }
}

library(reshape2)

reads2<-dcast(reads, Sample ~ Category, value.var = "Mapping_rate")
#head(reads2)
#subset <- reads2[-1]

#percentage_matrix <- apply(subset, 2, function(x) 100*x / subset$InputReads)

write.table(reads2, file = "quant_tests.tsv",quote=FALSE,sep="\t")
#write.table(percentage_matrix, file = "read_cleaning_results_percetage.txt",quote=FALSE,sep="\t")

