setwd("/mnt/disk5/2/3w3e/new_total/GPSeq/RadioHi-C/cal_equal_vol")

P_data <- data.frame()
for (i in 0:50) {
  name <- i
  P <- 1 - ((50-i)/50)^(1/3)
  aa <- c(i,P)
  P_data <- rbind(P_data,aa)
}
colnames(P_data) <- c("layer","P")
D2_score <- fread("/mnt/disk5/2/3w3e/new_total/GPSeq/RadioHi-C/Hs_100kb_2Dscore.bed",header = F) 
head(final_data)
final_data <- data.frame()
for (i in 1:(nrow(P_data)-1)) {
  layer_1 <- P_data[i+1,1]
  P_1 <- P_data[i,2]
  P_2 <- P_data[i+1,2]
  D2_score_tmp <- D2_score[D2_score$V4>=quantile(D2_score$V4,1-P_1) & D2_score$V4<=quantile(D2_score$V4,1-P_2),]
  D2_score_tmp$layer <- layer_1
  final_data <- rbind(final_data,D2_score_tmp)
}

write.table(final_data[,c(1,2,3,5)],file = "Hs_100kb_50layer_based_equal_vol.bed",col.names = F,row.names = F,quote = F,sep = "\t")

head(stats)
stats <- final_data %>% group_by(layer) %>% summarise(count=n())

ggplot(stats,aes(x=layer,y=count))+
  geom_bar(stat = "identity")+
  theme_bw()



