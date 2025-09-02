setwd("/mnt/disk6/3/useful_file/scHi-C/GAGE-seq")
library("data.table")
library("dplyr")

pairs_cell <- fread("cis_filter_K562_GSM7657692_contact_K3-0208-PL2-HiC.pairs",header = F)
bedfile <- read.table("100000_abs.bed",header = F)
chr_sd <- bedfile %>% group_by(V1) %>% summarise(min=min(V4),max=max(V4))

pairs_cell$V2 <- ceiling(pairs_cell$V2/100000)
pairs_cell$V4 <- ceiling(pairs_cell$V4/100000)

pairs_cell <- left_join(pairs_cell,chr_sd[,c(1,2)],by=c("V1"="V1"))
pairs_cell$V2 <- pairs_cell$V2 + pairs_cell$min - 1
pairs_cell$V4 <- pairs_cell$V4 + pairs_cell$min - 1
#pairs_cell[,c(2,4)] <- t(apply(pairs_cell[,c(2,4)],1,function(x) sort(x)))

cell_freq <- pairs_cell %>% group_by(V5,V2,V4) %>% summarise(freq=dplyr::n())
cell_freq <- cell_freq[cell_freq$V2 != cell_freq$V4,]
cell_freq_1 <- cell_freq[,c(1,3,2,4)]
colnames(cell_freq_1) <- colnames(cell_freq)
total <- rbind(cell_freq,cell_freq_1)
total$V5 <- gsub(",","_",total$V5)
#cell_newID <- as.data.frame(unique(sort(total$V5)))
write.table(total,file = "K562_GAGE_seq.freq_100kb.txt",col.names = F,row.names = F,sep = ",",quote = F)
#write.table(cell_newID,file = "cell_newID",col.names = F,row.names = F,quote = F)


