##call sub-boundary
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

bed_10kb <- read.table("/mnt/disk6/3/3w_manuscript/all_m3/10kb_matrix3/10000_abs.bed",header = F)

##以chr1为例
bed_10kb <- bed_10kb[bed_10kb$V1==args,]
bed_10kb$min <- bed_10kb$V4-15
bed_10kb$max <- bed_10kb$V4+15
bed_10kb <- bed_10kb[bed_10kb$min>0 & bed_10kb$max<max(bed_10kb$V4),]
bed_10kb$b_id <- paste("boundary",bed_10kb$V4,sep = "_")

output_file <- paste(args,"valid_bin.infomation",sep = "_")
write.table(bed_10kb,file = output_file,sep = "\t",col.names = T,row.names = F,quote = F)

##读取互作文件
matrix_10kb <- fread("/mnt/disk6/3/3w_manuscript/all_m3/10kb_matrix3/nor_Total_w3_28bp_M_10000.matrix3",header = F)
matrix_10kb_chr <- matrix_10kb[matrix_10kb$V1>=min(bed_10kb$V4) & matrix_10kb$V1<=max(bed_10kb$V4) &
                                 matrix_10kb$V2>=min(bed_10kb$V4) & matrix_10kb$V2<=max(bed_10kb$V4) & 
                                 matrix_10kb$V3>=min(bed_10kb$V4) & matrix_10kb$V3<=max(bed_10kb$V4),]
final_matrix <- data.frame()
for (i in 1:nrow(bed_10kb)) {
  aa <- bed_10kb$V4[i]
  matrix_10kb_1 <- matrix_10kb_chr
  matrix_10kb_1[,1:3] <- matrix_10kb_1[,1:3]-aa
  #matrix_10kb_1 <- matrix_10kb_1[abs(matrix_10kb_1$V1)<=15 & abs(matrix_10kb_1$V2)<=15 & abs(matrix_10kb_1$V3)<=15,]
  matrix_10kb_1 <- matrix_10kb_1[matrix_10kb_1$V1>= -15 & matrix_10kb_1$V1<=0 & 
                                   matrix_10kb_1$V2>=0 & matrix_10kb_1$V2<=15 & 
                                   matrix_10kb_1$V3>=0 & matrix_10kb_1$V3<=15,]
  matrix_10kb_1$b_id <- bed_10kb$b_id[i]
  final_matrix <- rbind(final_matrix,matrix_10kb_1)
}

IS_file <- paste(args,"3w_call_boundary_IC2.txt",sep = "_")
write.table(final_matrix,file = IS_file,col.names = F,row.names = F,quote = F,sep = "\t")

