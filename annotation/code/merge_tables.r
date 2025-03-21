library(dplyr, lib.loc="../data/")
sorghum <- read.table("../results/sorghum_blast.tbl", col.names = c("Sugarcane", "Arabidopsis", "evalue", "bitscore", "score"))[1:2]
maize <- read.table("../results/maize_blast.tbl", col.names = c("Sugarcane", "Maize", "evalue", "bitscore", "score"))[1:2] 
arabidopsis <- read.table("../results/arabidopsis_blast.tbl", col.names = c("Sugarcane", "Arabidopsis", "evalue", "bitscore", "score"))[1:2]

sorghum_arabidopsis <- full_join(sorghum, arabidopsis, by = "Sugarcane")
joined <- full_join(sorghum_arabidopsis, maize, by = "Sugarcane")
colnames(joined) <- c("Sugarcane", "Sorghum", "Arabidopsis", "Maize")

write.table(joined, "../results/sugarcane_sorghum_maize_arabidopsis.csv", quote = F, sep = ",", row.names = F, col.names = T)
