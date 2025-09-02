
imaging <- fread("/mnt/disk6/3/software/ChromatinImaging-master/Data/sort_probe_distances_parallel.csv2",header = F)

##top10
imaging_stats <- imaging %>%
  group_by(V2,V3) %>%
  slice_min(V4,n=10) %>%
  summarise(mean_top100 = mean(V4))


N_N_1 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/Cwalks/w2/new_shuf/chr21_N_N_1.txt",header = F)
N3_N3_2 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/Cwalks/w2/new_shuf/chr21_N3_N3_2.txt",header = F)
N4_N4_3 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/Cwalks/w2/new_shuf/chr21_N4_N4_3.txt",header = F)
Nk_Nk_k <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/Cwalks/w2/new_shuf/chr21_Nk_Nk_k.txt",header = F)

##Hi
N_N_1 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/HiPoreC_data/w2/HiPoreC_N_N_1.txt",header = F)
N3_N3_2 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/HiPoreC_data/w2/HiPoreC_N3_N3_2.txt",header = F)
N4_N4_3 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/HiPoreC_data/w2/HiPoreC_N4_N4_3.txt",header = F)
Nk_Nk_k <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/HiPoreC_data/w2/HiPoreC_Nk_Nk_k.txt",header = F)

N_N_1 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/mydata/2w/chr21_A_B.txt",header = F)
N3_N3_2 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/mydata/2w/chr21_A_G.txt",header = F)
N4_N4_3 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/mydata/2w/chr21_B_G.txt",header = F)

N_N_1 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/SPRITE_data/w2/chr21_SPRITE_A_B.txt",header = F)
N3_N3_2 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/SPRITE_data/w2/chr21_SPRITE_A_C.txt",header = F)
N4_N4_3 <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/SPRITE_data/w2/chr21_SPRITE_B_C.txt",header = F)


N_N_1[,1:2] <- N_N_1[,1:2] + 1
N_N_1[,1:2] <- data.frame(t(apply(N_N_1[,1:2], 1, function(x) sort(x))))
#N_N_1 <- N_N_1 %>% group_by(V1,V2) %>% summarise(count=dplyr::n())
N_N_1 <- inner_join(N_N_1,imaging_stats,by=c("V1"="V2","V2"="V3"))
#N_N_1_1 <- N_N_1[N_N_1$count>=quantile(N_N_1$count,0.98),]

N3_N3_2[,1:2] <- N3_N3_2[,1:2] + 1
N3_N3_2[,1:2] <- data.frame(t(apply(N3_N3_2[,1:2], 1, function(x) sort(x))))
#N3_N3_2 <- N3_N3_2 %>% group_by(V1,V2) %>% summarise(count=dplyr::n())
N3_N3_2 <- inner_join(N3_N3_2,imaging_stats,by=c("V1"="V2","V2"="V3"))

N4_N4_3[,1:2] <- N4_N4_3[,1:2] + 1
N4_N4_3[,1:2] <- data.frame(t(apply(N4_N4_3[,1:2], 1, function(x) sort(x))))
#N4_N4_3 <- N4_N4_3 %>% group_by(V1,V2) %>% summarise(count=dplyr::n())
N4_N4_3 <- inner_join(N4_N4_3,imaging_stats,by=c("V1"="V2","V2"="V3"))

Nk_Nk_k[,1:2] <- Nk_Nk_k[,1:2] + 1
Nk_Nk_k[,1:2] <- data.frame(t(apply(Nk_Nk_k[,1:2], 1, function(x) sort(x))))
#Nk_Nk_k <- Nk_Nk_k %>% group_by(V1,V2) %>% summarise(count=dplyr::n())
Nk_Nk_k <- inner_join(Nk_Nk_k,imaging_stats,by=c("V1"="V2","V2"="V3"))

Cwalks_total <- as.data.frame(c(N_N_1$mean_top100,N3_N3_2$mean_top100,N4_N4_3$mean_top100,Nk_Nk_k$mean_top100))
Cwalks_total$Type <- c(rep("N1",nrow(N_N_1)),rep("N2",nrow(N3_N3_2)),rep("N3",nrow(N4_N4_3)),rep("N4",nrow(Nk_Nk_k)))
colnames(Cwalks_total) <- c("dis","Type")

Hi_total <- as.data.frame(c(N_N_1$mean_top100,N3_N3_2$mean_top100,N4_N4_3$mean_top100,Nk_Nk_k$mean_top100))
Hi_total$Type <- c(rep("N1",nrow(N_N_1)),rep("N2",nrow(N3_N3_2)),rep("N3",nrow(N4_N4_3)),rep("N4",nrow(Nk_Nk_k)))
colnames(Hi_total) <- c("dis","Type")

sprite_total <- as.data.frame(c(N_N_1$mean_top100,N3_N3_2$mean_top100,N4_N4_3$mean_top100))
sprite_total$Type <- c(rep("N1",nrow(N_N_1)),rep("N2",nrow(N3_N3_2)),rep("N3",nrow(N4_N4_3)))
colnames(sprite_total) <- c("dis","Type")

w3_total <- as.data.frame(c(N_N_1$mean_top100,N3_N3_2$mean_top100,N4_N4_3$mean_top100))
w3_total$Type <- c(rep("N1",nrow(N_N_1)),rep("N2",nrow(N3_N3_2)),rep("N3",nrow(N4_N4_3)))
colnames(w3_total) <- c("dis","Type")

write.table(Total,file = "/mnt/disk6/3/3w_manuscript/fig2/FISH/Total_2w_dis_top10_cell_all.txt",quote = F,col.names = F,row.names = F,sep = "\t")

Total <- rbind(Cwalks_total,Hi_total,sprite_total,w3_total)
Total$Method <- c(rep("C-walks",nrow(Cwalks_total)),
                  rep("HiPore-C",nrow(Hi_total)),
                  rep("SPRITE",nrow(sprite_total)),
                  rep("3wHi-C",nrow(w3_total)))

Total <- read.table("/mnt/disk6/3/3w_manuscript/fig2/FISH/Total_2w_dis_top10_cell_all.txt",header = F)
colnames(Total) <- c("dis","Type","Method")
mean <- Total %>% group_by(Type,Method) %>% summarise(mean=mean(dis))
sd <- Total %>% group_by(Type,Method) %>% summarise(sd=sd(dis))
df_res <- data.frame(mean, sd=sd$sd)
colnames(df_res) = c("Species","Method","Mean","Sd")
write.table(df_res,file = "/mnt/disk6/3/3w_manuscript/fig2/FISH/top10_dis_diff_methods.txt",col.names = T,row.names = F,quote = F,sep = "\t")

p1 <- ggplot(Total,aes(x=Type,y=dis,fill=Method))+
  theme_bw()+
  geom_violin(position=position_dodge(width=0.8),size=0.75)+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),size=0.3,outlier.colour = NA)+
  #geom_bar(stat="identity", position=position_dodge(),color="black", width=.5)+
  #geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd),position=position_dodge(0.5), width=.2,color="#3B3838")+
  scale_fill_manual(values = c('#427ab2','#48c0aa','#f6c957','#AD8BC6'))+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        #axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black", size = 0.5),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.text.y=element_text(family="Times",size=17,colour="black"),
        axis.text.x=element_text(family="Times",size=17,colour="black"),
        axis.title.x=element_text(family="Times",size = 17,colour="black"),
        axis.title.y=element_text(family="Times",size = 17,colour="black"),
        legend.position = "top",legend.text = element_text(family="Times",size = 15,colour="black"),
        legend.title = element_blank())+
  xlab("")+ylab("Mean spatial distance(nm)")
ggsave("/mnt/disk6/3/software/ChromatinImaging-master/Data/mydata/2w/2w_dis.pdf", plot = p1, dpi = 300,width = 3.1,height = 3)

df_res <- read.table("/mnt/disk6/3/3w_manuscript/fig2/FISH/top10_dis_diff_methods.txt",header = T)


p1 <- ggplot(df_res,aes(x=Method,y=Mean,fill=Species))+
  theme_bw()+
  #geom_boxplot()+
  geom_bar(stat="identity", position=position_dodge(0.8),color="black", width=.6)+
  #geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd),position=position_dodge(0.8), width=.2,color="#3B3838")+
  scale_fill_manual(values = c('#386CAF','#7FC87F','#fc9533','#9f79d3'))+ #'#427ab2','#48c0aa','#f6c957','#AD8BC6'
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        #axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black", size = 0.5),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.text.y=element_text(family="Times",size=17,colour="black"),
        axis.text.x=element_text(family="Times",size=17,colour="black"),
        axis.title.x=element_text(family="Times",size = 17,colour="black"),
        axis.title.y=element_text(family="Times",size = 17,colour="black"),
        legend.position = "none",
        #legend.position = "top",legend.text = element_text(family="Times",size = 15,colour="black"),
        legend.title = element_blank())+
  xlab("")+ylab("Spatial distance (nm)")+ylim(0,22)   #ylim(0,600)
p1

ggsave("/mnt/disk6/3/3w_manuscript/fig2/FISH/ALL_dis_top10_mean.pdf", plot = p1, dpi = 300,width = 5,height = 2.8)











