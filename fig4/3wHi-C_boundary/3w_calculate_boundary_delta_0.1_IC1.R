##call 3w boundary
library(data.table)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
boundary_file <- paste(args,"3w_call_boundary_IC1.txt",sep = "_")
TAD_boundary_all <- fread(boundary_file,header = F)

TAD_boundary_all <- TAD_boundary_all[TAD_boundary_all$V1 != 0 & TAD_boundary_all$V2 != 0 & TAD_boundary_all$V3 != 0,]
TAD_boundary_all_1 <- TAD_boundary_all[,c(1,3,2,4,5)]
TAD_boundary_all_2 <- TAD_boundary_all[,c(2,1,3,4,5)]
TAD_boundary_all_3 <- TAD_boundary_all[,c(2,3,1,4,5)]
TAD_boundary_all_4 <- TAD_boundary_all[,c(3,1,2,4,5)]
TAD_boundary_all_5 <- TAD_boundary_all[,c(3,2,1,4,5)]

colnames(TAD_boundary_all_1) <- colnames(TAD_boundary_all)
colnames(TAD_boundary_all_2) <- colnames(TAD_boundary_all)
colnames(TAD_boundary_all_3) <- colnames(TAD_boundary_all)
colnames(TAD_boundary_all_4) <- colnames(TAD_boundary_all)
colnames(TAD_boundary_all_5) <- colnames(TAD_boundary_all)

TAD_boundary_all <- rbind(TAD_boundary_all,TAD_boundary_all_1,TAD_boundary_all_2,
                          TAD_boundary_all_3,TAD_boundary_all_4,
                          TAD_boundary_all_5)
rm(TAD_boundary_all_1,TAD_boundary_all_2,
   TAD_boundary_all_3,TAD_boundary_all_4,
   TAD_boundary_all_5)

TAD_boundary_all_clean <- TAD_boundary_all[!(TAD_boundary_all$V1 == TAD_boundary_all$V2 &
                                               TAD_boundary_all$V1 == TAD_boundary_all$V3 &
                                               TAD_boundary_all$V2 == TAD_boundary_all$V3),]

boundary_statss <- TAD_boundary_all_clean %>% group_by(V5) %>% summarise(freq=sum(V4))
boundary_statss$b_num <- as.numeric(substring(boundary_statss$V5,10))


valid_file <- paste(args,"valid_bin.infomation",sep = "_")
valid_b <- fread(valid_file,header = T)

#valid_b <- read.table("chr18_valid_bin.infomation",header = T)
boundary_statss$signal <- log2(boundary_statss$freq/mean(boundary_statss$freq))

valid_b_1 <- inner_join(valid_b[,c(1:4,7)],boundary_statss[,c(1,4)],by=c("b_id"="V5"))

#10
#4
valid_b_1$left <- sapply(valid_b_1$V4, function(x) {
  left_rows <- valid_b_1[valid_b_1$V4>=(x-10) & valid_b_1$V4<x,]
  sum(left_rows$signal)/10
})

valid_b_1$right <- sapply(valid_b_1$V4, function(x) {
  right_rows <- valid_b_1[valid_b_1$V4>x & valid_b_1$V4<=(x+10),]
  sum(right_rows$signal)/10
})

valid_b_1$delta <- valid_b_1$left - valid_b_1$right


final_valid_1 <- data.frame()
for (i in 1:(nrow(valid_b_1)-1)) {
  if (valid_b_1[i,4]>=(min(valid_b_1$V4)+10) & valid_b_1[i,4]<=(max(valid_b_1$V4)-10) & 
      valid_b_1[i,9]*valid_b_1[i+1,9]<0 & valid_b_1[i,9]>0 & valid_b_1[i+1,9]<0) {
    new_row <- valid_b_1[i,]
    final_valid_1 <- rbind(final_valid_1,new_row)
  }
}

final_valid_2 <- data.frame()
for (i in 1:(nrow(valid_b_1)-1)) {
  subset_left <- valid_b_1[valid_b_1$V4 >= valid_b_1[i,]$V4 - 10 & valid_b_1$V4 <= valid_b_1[i,]$V4,9]
  subset_right <- valid_b_1[valid_b_1$V4 >= valid_b_1[i,]$V4 & valid_b_1$V4 <= valid_b_1[i,]$V4 + 10,9]
  if (valid_b_1[i,4]>=(min(valid_b_1$V4)+10) & valid_b_1[i,4]<=(max(valid_b_1$V4)-10) & 
      valid_b_1[i,9]*valid_b_1[i+1,9]<0 & valid_b_1[i,9]>0 & valid_b_1[i+1,9]<0 & valid_b_1[i+1,4]-valid_b_1[i,4]==1 & 
      max(subset_left) - min(subset_right) > 0.1) {
    new_row <- valid_b_1[i,]
    final_valid_2 <- rbind(final_valid_2,new_row)
  }
}
IS_file <- paste("/mnt/disk6/3/3w_manuscript/sub_boundary/method2_call_boundary/ice_is-300Kb-ids-200kb_right_up/",args,"_3w_IS_RU.txt",sep = "")
delta_file <- paste("/mnt/disk6/3/3w_manuscript/sub_boundary/method2_call_boundary/ice_is-300Kb-ids-200kb_right_up/",args,"_3w_delta_RU_gt0.1.txt",sep = "")
final_boundary_file <- paste("/mnt/disk6/3/3w_manuscript/sub_boundary/method2_call_boundary/ice_is-300Kb-ids-200kb_right_up/",args,"_3w_true_boundary_RU_gt0.1.txt",sep = "")
final_boundary_file2 <- paste("/mnt/disk6/3/3w_manuscript/sub_boundary/method2_call_boundary/ice_is-300Kb-ids-200kb_right_up/",args,"_3w_true_boundary_RU.txt",sep = "")

write.table(valid_b_1[,c(1:4,6)],file = IS_file,col.names = T,row.names = F,quote = F,sep = "\t")
write.table(valid_b_1[,c(1:4,9)],file = delta_file,col.names = F,row.names = F,quote = F,sep = "\t")
write.table(final_valid_2[,1:4],file = final_boundary_file,col.names = F,row.names = F,quote = F,sep = "\t")
write.table(final_valid_1[,1:4],file = final_boundary_file2,col.names = F,row.names = F,quote = F,sep = "\t")










