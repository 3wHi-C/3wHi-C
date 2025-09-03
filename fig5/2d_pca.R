#!/usr/bin/env Rscript
library(dplyr)
library(doParallel)
library(foreach)


#========PATH==========
matrix_500Kb <- read.table("/mnt/disk6/3/3w_manuscript/all_m3/250kb_matrix3/OE/chr16_OE_250kb.matrix3",header = F,sep = "\t")
##========input AB infomation=========
#/mnt/disk5/2/3w3e/new_total/comp_bin/250Kb/PC1_AB_noNA.txt
#/mnt/disk5/2/3w3e/new_total/comp_bin/500Kb/PC1_AB_noNA.txt
ABcomp <- read.table("/mnt/disk5/2/3w3e/new_total/comp_bin/250Kb/PC1_AB_noNA.txt",header = F)

ABcomp$binID <- as.numeric(substring(ABcomp$V5,3))

matrix_500Kb <- inner_join(matrix_500Kb,ABcomp[,c(6,5)],by=c("V1"="binID"))
matrix_500Kb <- inner_join(matrix_500Kb,ABcomp[,c(6,5)],by=c("V2"="binID"))
matrix_500Kb <- inner_join(matrix_500Kb,ABcomp[,c(6,5)],by=c("V3"="binID"))

colnames(matrix_500Kb)[5:7] <- c("bin1","bin2","bin3")

##==========Completion matrix==========
#matrix_250Kb$length <- matrix_250Kb$V3 - matrix_250Kb$V1

matrix_500Kb_1 <- matrix_500Kb[,c(1,3,2,4,5,7,6)]
matrix_500Kb_2 <- matrix_500Kb[,c(2,1,3,4,6,5,7)]
matrix_500Kb_3 <- matrix_500Kb[,c(2,3,1,4,6,7,5)]
matrix_500Kb_4 <- matrix_500Kb[,c(3,1,2,4,7,5,6)]
matrix_500Kb_5 <- matrix_500Kb[,c(3,2,1,4,7,6,5)]
colnames(matrix_500Kb_1) <- colnames(matrix_500Kb)
colnames(matrix_500Kb_2) <- colnames(matrix_500Kb)
colnames(matrix_500Kb_3) <- colnames(matrix_500Kb)
colnames(matrix_500Kb_4) <- colnames(matrix_500Kb)
colnames(matrix_500Kb_5) <- colnames(matrix_500Kb)

total_matrix_500Kb <- rbind(matrix_500Kb,matrix_500Kb_1,matrix_500Kb_2,matrix_500Kb_3,matrix_500Kb_4,matrix_500Kb_5)
rm(matrix_500Kb_1,matrix_500Kb_2,matrix_500Kb_3,matrix_500Kb_4,matrix_500Kb_5)

##========AA | BB================
comp <- c("A","B")
AA_BB <- data.frame()
for (i in comp) {
  matrix_500Kb_2 <- total_matrix_500Kb %>% filter(grepl(i,bin1)) %>% filter(grepl(i,bin2))
  AA_BB <- rbind(AA_BB,matrix_500Kb_2)
}

AA_BB[,5:6] <- data.frame(t(apply(AA_BB[,5:6],1,function(x) sort(x))))
AA_BB$bin12 <- paste(AA_BB$bin1,AA_BB$bin2,sep = "-")
total_bin <- as.data.frame(c(unique(sort(AA_BB$bin12))))
total_bin$ID <- 1:nrow(total_bin)
colnames(total_bin) <- c("pairs","ID")
write.table(total_bin,file = "/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/pairs_ID.txt",sep="\t",row.names = F,col.names = F,quote = F)

library("data.table")
cor_total <- fread("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/Total_chr16_signaficant.txt",header = F)
colnames(cor_total) <- c("ID1","ID2","bin1","bin2","cor","P")
N <- nrow(total_bin)

cor_total <- cor_total[cor_total$bin1 %in% total_bin$pairs & cor_total$bin2 %in% total_bin$pairs,]

cor_total <- cor_total %>% left_join(total_bin,by=c("bin1"="pairs")) %>% left_join(total_bin,by=c("bin2"="pairs")) 
write.table(cor_total[,c(7,8,5)],file = "/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/chr16_sparse_matrix.txt",sep="\t",row.names = F,col.names = F,quote = F)

cor_total <- fread("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/Total_chr16_signaficant.txt",header = F)
colnames(cor_total) <- c("ID1","ID2","bin1","bin2","cor","P")
N <- nrow(total_bin)

oe_matrix_cor <- matrix(0, nrow = N, ncol = N,dimnames = list(total_bin$ID,total_bin$ID))

cor_total2 <- cor_total[,c(2,1,4,3,5,6)]
colnames(cor_total2) <- colnames(cor_total)
total_cor_chr20_matrix <- rbind(cor_total,cor_total2)
rm(cor_total2)
total_cor_chr20_matrix <- total_cor_chr20_matrix %>% left_join(total_bin,by=c("bin1"="pairs")) %>% left_join(total_bin,by=c("bin2"="pairs")) 

for (t in 1:nrow(total_cor_chr20_matrix)) {
  oe_matrix_cor[total_cor_chr20_matrix[t,]$ID.x,total_cor_chr20_matrix[t,]$ID.y] <- total_cor_chr20_matrix[t,]$cor
}

write.table(oe_matrix_cor,file="/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/chr16_matrix.txt",sep="\t",row.names = T,col.names = T,quote = F)

##=================================first-step dimension reduction========================================
oe_matrix_cor <- read.table("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/chr16_matrix.txt",header = T,row.names=1,sep = "\t")
colnames(oe_matrix_cor) <- 1:ncol(oe_matrix_cor)
#oe_matrix_cor <- matrix(oe_matrix_cor)
pca_result <- prcomp(oe_matrix_cor)
PC_values <- as.data.frame(pca_result$x[,1:10])
chromosome <- ABcomp[ABcomp$V1=="chr16",]
length_m <- nrow(chromosome)
matrix_name <- c(chromosome$V5)
PC_values$pairs <- total_bin$pairs
for (i in 1:nrow(PC_values)) {
  PC_values[i,12] <- strsplit(PC_values[i,11],split="-")[[1]][1]
  PC_values[i,13] <- strsplit(PC_values[i,11],split="-")[[1]][2]
}

write.table(PC_values,file = "/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/chr16_PC1_value_pairs.txt",col.names = T,row.names = T,sep = "\t",quote = F)


##===============================  second-step dimension reduction=======================================
tianchong_matrix <- PC_values[,c(12,13,1)]
colnames(tianchong_matrix) <- c("bin1","bin2","PC1")
tianchong_matrix2 <- PC_values[,c(13,12,1)]
colnames(tianchong_matrix2) <- c("bin1","bin2","PC1")
total_tianchong <- rbind(tianchong_matrix,tianchong_matrix2)
rm(tianchong_matrix2)

PC1_matrix <- matrix(0, nrow = length_m, ncol = length_m,dimnames = list(matrix_name,matrix_name))
for (t in 1:nrow(total_tianchong)) {
  PC1_matrix[total_tianchong[t,]$bin1,total_tianchong[t,]$bin2] <- total_tianchong[t,]$PC1
}
write.table(PC1_matrix,file="/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/chr16_1D.matrix.txt",col.names = F,row.names = F,quote=F,sep = "\t")

pc1_pca <- prcomp(PC1_matrix,center=T,scale=T)
pc1_pca_values <- pc1_pca$x[,1:10]
pc1_pca_values <- as.data.frame(pc1_pca_values)
pc1_pca_values$bin <- rownames(pc1_pca_values)
pc1_pca_values <- left_join(pc1_pca_values,ABcomp[,4:5],by=c("bin"="V5"))
cor.test(pc1_pca_values$PC1,pc1_pca_values$V4)

write.table(pc1_pca_values[,c(11,1)],file = "/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/chr16_1D_PC1.txt",sep = "\t",col.names = F,row.names = F,quote = F)









