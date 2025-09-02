##======================================figS2c: distance==============================================
poreC <- read.table("/mnt/disk1/1/gy/HiPore-C/distance/pic/Pore-C-chr.1.distance.txt",header = F)
poreC <- poreC[rowSums(poreC == 0) == 0, ]

cwalks <- read.table("/mnt/disk1/1/gy/Cwalks/distance/intra_chr1.Cwalks.distance.txt",header = F)
cwalks <- cwalks[rowSums(cwalks == 0) == 0, ]

three_way <- fread("/mnt/disk1/1/gy/HiPore-C/distance/pic/three-way-01-02-12.txt",header = F)
head(three_way)
three_way <- three_way[rowSums(three_way == 0) == 0, ]

sprite_1 <- fread("/mnt/disk1/1/gy/SPRITE/GSE114242_human_combined_clusters/hg38_final_3w/final_hg38_file1.bed",header = F)
sprite_2 <- fread("/mnt/disk1/1/gy/SPRITE/GSE114242_human_combined_clusters/hg38_final_3w/final_hg38_file2.bed",header = F)
sprite_3 <- fread("/mnt/disk1/1/gy/SPRITE/GSE114242_human_combined_clusters/hg38_final_3w/final_hg38_file3.bed",header = F)
sprite_final <- sprite_1 %>% inner_join(sprite_2,by=c("V4"="V4")) %>% inner_join(sprite_3,by=c("V4"="V4"))
sprite_final <- sprite_final[sprite_final$V1.x==sprite_final$V1.y & sprite_final$V1.x==sprite_final$V1 & 
                               sprite_final$V1.y==sprite_final$V1 & sprite_final$V1.x=="chr1",]

sprite_final$AB <- abs(sprite_final$V2.x - sprite_final$V2.y)
sprite_final$AC <- abs(sprite_final$V2.x - sprite_final$V2)
sprite_final$BC <- abs(sprite_final$V2.y - sprite_final$V2)

CHIA_Drop <- read.table("/mnt/disk6/3/useful_file/CHIA-Drop_Drosophila/shuf_pos_w3.bed",header = F)
CHIA_Drop$AB <- abs(CHIA_Drop$V2 - CHIA_Drop$V3)
CHIA_Drop$AC <- abs(CHIA_Drop$V2 - CHIA_Drop$V4)
CHIA_Drop$BC <- abs(CHIA_Drop$V3 - CHIA_Drop$V4)

##Tri-C
Tric <- fread("/mnt/disk1/1/gy/Tri-C/dis_C57_GSE107755_b_CTCF_Tri-C_Samples_Combined.txt",header = F)
head(Tric)
Tric$dis <- abs(Tric$V2-Tric$V3)
merged <- c(abs(poreC$V1),abs(poreC$V2),abs(poreC$V3),
            abs(cwalks$V1),abs(cwalks$V2),abs(cwalks$V3),
            abs(three_way$V1),abs(three_way$V3),abs(three_way$V2),
            abs(sprite_final$AB),abs(sprite_final$AC),abs(sprite_final$BC),CHIA_Drop$AB,CHIA_Drop$AC,CHIA_Drop$BC,Tric$dis)

merged_final <- data.frame(MergedColumn = merged, DataFrame = rep(c("HiPore-C","HiPore-C","HiPore-C",
                                                                    "C-walks","C-walks","C-walks",
                                                                    "3wHi-C","3wHi-C","3wHi-C",
                                                                    "SPRITE","SPRITE","SPRITE",
                                                                    "ChIA-Drop","ChIA-Drop","ChIA-Drop","Tri-C"),
                                                                  times = c(nrow(poreC),nrow(poreC),nrow(poreC),
                                                                            nrow(cwalks),nrow(cwalks),nrow(cwalks),
                                                                            nrow(three_way),nrow(three_way),nrow(three_way),
                                                                            nrow(sprite_final),nrow(sprite_final),nrow(sprite_final),
                                                                            nrow(CHIA_Drop),nrow(CHIA_Drop),nrow(CHIA_Drop),nrow(Tric))))

colnames(merged_final) <- c("distance","Type")

library("RColorBrewer")
cl <- brewer.pal(9,"Set1")

#merged_final_1 <- merged_final[c(1:1000,400000:4000,1100000:1100100),]
library("ggplot2")
pdf("/mnt/disk6/3/3w_manuscript/fig2/Cwalks_dis.pdf", width = 8, height = 6.5, bg = "transparent")
p1 <- ggplot(merged_final,aes(x = log10(distance),color=Type))+
  stat_ecdf(geom="smooth",se=F,linewidth=1.5)+
  scale_color_manual(values = cl)+
  ylab("Cumulative frequency")+
  xlab("log10(distance(bp))")+
  theme_bw()+
  geom_hline(yintercept = 0.5,linetype = "dashed", color = "black")+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.text.y=element_text(family="Times",size=17,colour="black"),
        axis.text.x=element_text(family="Times",size=17,colour="black"),
        axis.title.x=element_text(family="Times",size = 17,colour="black"),
        axis.title.y=element_text(family="Times",size = 17,colour="black"),
        legend.text=element_text(family="Times", colour="black",size=15),
        legend.title = element_blank())+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.05))+
  scale_x_continuous(expand = c(0, 0),limits = c(0,8.5),breaks = c(2,4,6,8))

#ggsave("/mnt/disk6/3/3w_manuscript/fig2/distance/SPRITE_all_dis.pdf", plot = p1, dpi = 300,width = 8,height = 6.5)
ggsave("/mnt/disk6/3/3w_manuscript/fig2/distance/dis_6_method_All.tiff", plot = p1, dpi = 300,width = 11,height = 6)

##================================fig2c: Radial-C==============================================================
#/mnt/disk5/2/3w3e/new_total/GPSeq/HiPoreC/HiPoreC_N3_N3+2_layer_diff.txt
poreC_1 <- read.table("/mnt/disk5/2/3w3e/new_total/GPSeq/HiPoreC/50_layers_equal_vol/HiPoreC_N3_N3+2_layer_diff.txt",header = F)
poreC_2 <- read.table("/mnt/disk5/2/3w3e/new_total/GPSeq/HiPoreC/50_layers_equal_vol/HiPoreC_N4_N4+3_layer_diff.txt",header = F)
poreC_3 <- read.table("/mnt/disk5/2/3w3e/new_total/GPSeq/HiPoreC/50_layers_equal_vol/HiPoreC_NK_NK+K_layer_diff.txt",header = F)
poreC_4 <- read.table("/mnt/disk5/2/3w3e/new_total/GPSeq/HiPoreC/50_layers_equal_vol/HiPoreC_N_N+1_layer_diff.txt",header = F)

C_1 <- read.table("/mnt/disk5/2/3w3e/new_total/GPSeq/C-walks/50_layers_equal_vol/new_shuf/Cwalks_N_N+1_layer_diff.txt",header = F)
C_2 <- read.table("/mnt/disk5/2/3w3e/new_total/GPSeq/C-walks/50_layers_equal_vol/new_shuf/Cwalks_N3_N3+2_layer_diff.txt",header = F)
C_3 <- read.table("/mnt/disk5/2/3w3e/new_total/GPSeq/C-walks/50_layers_equal_vol/new_shuf/Cwalks_N4_N4+3_layer_diff.txt",header = F)
C_4 <- read.table("/mnt/disk5/2/3w3e/new_total/GPSeq/C-walks/50_layers_equal_vol/new_shuf/Cwalks_Nk_Nk+k_layer_diff.txt",header = F)

w3 <- fread("/mnt/disk5/2/3w3e/new_total/GPSeq/3wHi-C/50_layers_equal_vol/3wHi-C_layer_diff.txt",header = F)
sprite <- fread("/mnt/disk5/2/3w3e/new_total/GPSeq/SPRITE/50_layers_equal_vol/SPRITE_layer_diff.txt",header = F)
head(sprite)

merged <- c(
            #abs(w3$V1),abs(w3$V2),abs(w3$V3))
            #abs(poreC_4$V1),abs(poreC_1$V1),abs(poreC_2$V1),abs(poreC_3$V1))
            #abs(sprite$V1),abs(sprite$V2),abs(sprite$V3))
            abs(C_1$V1),abs(C_2$V1),abs(C_3$V1),abs(C_4$V1))
#"3wHi-C-Pairs","3wHi-C-Pairs","3wHi-C-Pairs"
merged_final <- data.frame(MergedColumn = merged,
                           DataFrame = rep(c(
                                            #"3wHi-C-α/β","3wHi-C-α/γ","3wHi-C-β/γ"),
                                             #"HiPoreC-N/N+1","HiPoreC-N/N+2","HiPoreC-N/N+3","HiPoreC-N/N+k"),
                                             #"A/B","A/C","B/C"),
                                            "C-walks-N/N+1","C-walks-N/N+2","C-walks-N/N+3","C-walks-N/N+k"),
                                           times = c(
                                             #nrow(w3),nrow(w3),nrow(w3))))
                                                     #nrow(poreC_4),nrow(poreC_1),nrow(poreC_2),nrow(poreC_3))))
                                                     #nrow(sprite),nrow(sprite),nrow(sprite))))
                                             nrow(C_1),nrow(C_2),nrow(C_3),nrow(C_4))))

colnames(merged_final) <- c("Layer_Difference","Type")

library("RColorBrewer")
library("ggsci")
coul <- colorRampPalette(brewer.pal(9,"Paired"))(9)
cl <- pal_lancet("lanonc",alpha = 0.9)(6)
p2 <- ggplot(merged_final,aes(x = Layer_Difference,color=Type))+ #,linetype=Type
  stat_ecdf(geom="smooth",se=F,linewidth=1.5)+
  scale_color_manual(values = c('#386CAF','#7FC87F','#fc9533','#9f79d3'))+
  scale_linetype_manual(values = c("solid","longdash","dotted","twodash"))+
  ylab("Cumulative differences of layers %")+
  xlab("Layer difference (cis) in different methods")+
  theme_bw()+
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
        legend.text=element_text(family="Times", colour="black",size=17),
        legend.title = element_blank(),legend.position = c(0.85,0.1),legend.justification = c(0.85,0.1))+ #legend.position = c(0.85,0.1)
  scale_y_continuous(expand = c(0,0),limits = c(0,1.05))+
  scale_x_continuous(limits = c(0,50))


ggsave("/mnt/disk6/3/3w_manuscript/fig2/50_layers_equal_vol/C-walks_layer_tmp.pdf", plot = p2, dpi = 300,width = 5,height = 5)

##==================================fig2d: compare 3wHi-C and SPRITE=======================================
w3 <- fread("/mnt/disk5/2/3w3e/new_total/GPSeq/3wHi-C/50_layers_equal_vol/3wHi-C_layer_diff.txt",header = F)
sprite <- fread("/mnt/disk5/2/3w3e/new_total/GPSeq/SPRITE/50_layers_equal_vol/SPRITE_layer_diff.txt",header = F)
merged <- c(
  abs(w3$V1),abs(w3$V2),abs(w3$V3),abs(sprite$V1),abs(sprite$V2),abs(sprite$V3))
merged_final <- data.frame(MergedColumn = merged,
                           DataFrame = rep(c("3wHi-C","3wHi-C","3wHi-C","SPRITE","SPRITE","SPRITE"),
                             times = c(nrow(w3),nrow(w3),nrow(w3),nrow(sprite),nrow(sprite),nrow(sprite))))
head(merged_final)
merged_final$Type <- c(rep("A/B",nrow(w3)),rep("A/C",nrow(w3)),rep("B/C",nrow(w3)),
                       rep("A/B",nrow(sprite)),rep("A/C",nrow(sprite)),rep("B/C",nrow(sprite)))
#colnames(merged_final) <- c("layer","method","type")
colnames(merged_final) <- c("layer","method")
merged_final$type2 <- 
  ifelse(merged_final$layer==0,"0",ifelse(merged_final$layer>=1 & merged_final$layer<=2,"[1,2]",
                                          ifelse(merged_final$layer>=3 & merged_final$layer<=5,"[3,5]",
                                  ifelse(merged_final$layer>=6 & merged_final$layer<=10,"[6,10]",
                                         ifelse(merged_final$layer>=11 & merged_final$layer<=20,"[11,20]","[21,50]")))))

stats <- merged_final %>% group_by(method,type2) %>% summarise(count=n())
stats_2 <- stats %>% group_by(method) %>% mutate(percent=count/sum(count))

aa <- factor(stats_2$type2,levels = c("0","[1,2]","[3,5]","[6,10]","[11,20]","[21,50]"))
p1 <- ggplot(data=stats_2,aes(x=aa,y=percent,fill=method))+
  geom_bar(stat="identity",position=position_dodge(0.8),width=0.8)+
  scale_fill_manual(values = c('#427ab2','#f6c957','#48c0aa')) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, colour = "black", size = 15, angle = 0),
        axis.text.y = element_text(size = 16, face = "plain", colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 2, hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        #axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black", size = 0.5),
        axis.line.x = element_line(colour = "black", size = 0.5),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "lines"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "plain", colour = "black"),
        plot.title = element_text(colour = "black", face = "plain", size = 15, vjust = 8, hjust = 1))+
  ylim(0,0.5)

ggsave("/mnt/disk6/3/3w_manuscript/fig2/50_layers_equal_vol/3wHi-C_and_sprite_layer_bar.pdf",
       plot = p1, dpi = 300,width = 4,height = 2.9)

##========================================figS2d: 1D layer==========================================   

w3 <- fread("/mnt/disk5/2/3w3e/new_total/GPSeq/3wHi-C/50_layers_equal_vol/1D/1D_layer",header = F)
hi <- fread("/mnt/disk5/2/3w3e/new_total/GPSeq/HiPoreC/50_layers_equal_vol/1D/1D_layer",header = F)
sprite <- fread("/mnt/disk5/2/3w3e/new_total/GPSeq/SPRITE/50_layers_equal_vol/1D/1D_layer",header = F)
cwalks <- fread("/mnt/disk5/2/3w3e/new_total/GPSeq/C-walks/50_layers_equal_vol/1D/1D_layer",header = F)

total <- rbind(w3,hi,sprite,cwalks)
total$method <- c(rep("3wHi-C",nrow(w3)),rep("HiPore-C",nrow(hi)),rep("SPRITE",nrow(sprite)),rep("C-walks",nrow(cwalks)))
total$type2 <- ifelse(total$V1>=1 & total$V1<=2,"[1,2]",
                      ifelse(total$V1>=3 & total$V1<=5,"[3,5]",
                             ifelse(total$V1>=6 & total$V1<=10,"[6,10]",
                                    ifelse(total$V1>=11 & total$V1<=20,"[11,20]","[21,50]"))))

total_1 <- total %>% group_by(method,type2) %>% summarise(count=n())
total_1 <- total_1 %>% group_by(method) %>% mutate(percent=count/sum(count))

aa <- factor(total_1$type2,levels = c("[1,2]","[3,5]","[6,10]","[11,20]","[21,50]"))
  p1 <- ggplot(data=total_1,aes(x=aa,y=percent,fill=method))+
    geom_bar(stat="identity",position=position_dodge(0.8),width=0.8)+
    scale_fill_manual(values = c('#386CAF','#7FC87F','#fc9533','#9f79d3')) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, colour = "black", size = 15, angle = 0),
          axis.text.y = element_text(size = 16, face = "plain", colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 2, hjust = 0.5),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.y = element_line(colour = "black", size = 0.5),
          axis.ticks.x = element_line(colour = "black", size = 0.5),
          #axis.ticks.x = element_blank(),
          axis.line.y = element_line(colour = "black", size = 0.5),
          axis.line.x = element_line(colour = "black", size = 0.5),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "lines"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 12, face = "plain", colour = "black"),
          plot.title = element_text(colour = "black", face = "plain", size = 15, vjust = 8, hjust = 1))+
    ylim(0,0.91)
  
ggsave("/mnt/disk6/3/3w_manuscript/fig2/50_layers_equal_vol/1D_layer.pdf", plot = p1, dpi = 300,width = 4,height = 3) 

##============================================fig2b===============================================================
H1 <- fread("/mnt/disk1/1/gy/HiPore-C/distance/chr1-N2_N2+1.txt",header = F)
H2 <- fread("/mnt/disk1/1/gy/HiPore-C/distance/chr1-N_N+2.txt",header = F)
H3 <- fread("/mnt/disk1/1/gy/HiPore-C/distance/chr-N4_N4+3.txt",header = F)
H4 <- fread("/mnt/disk1/1/gy/HiPore-C/distance/chr1-Nk_Nk+k.txt",header = F)

H1 <- H1 %>% sample_n(100000)
H2 <- H2 %>% sample_n(100000)
H3 <- H3 %>% sample_n(100000)
H4 <- H4 %>% sample_n(100000)
merged <- c(abs(H1$V2),abs(H2$V2),abs(H3$V2),abs(H4$V2))
merged_final <- data.frame(MergedColumn = merged, DataFrame = rep(c("N/N+1","N/N+2","N/N+3","N/N+k"),
                                                                  times = c(nrow(H1),nrow(H2),nrow(H3),nrow(H4))))
head(merged_final)
colnames(merged_final) <- c("dis","type")
p2 <- ggplot(merged_final,aes(x = log10(dis),color=type))+
  stat_ecdf(geom="smooth",se=F,linewidth=3)+
  scale_color_manual(values = c('#386CAF','#7FC87F','#fc9533','#9f79d3'))+ 
  ylab("Cumulative frequency")+
  xlab("log10(distance(bp))")+
  theme_bw()+
  geom_hline(yintercept = 0.5,linetype = "dashed", color = "black")+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.ticks.x = element_line(colour = "black", size = 1),
        #axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.text.y=element_text(family="Times",size=17,colour="black"),
        axis.text.x=element_text(family="Times",size=17,colour="black"),
        axis.title.x=element_text(family="Times",size = 17,colour="black"),
        axis.title.y=element_text(family="Times",size = 17,colour="black"),
        legend.text=element_text(family="Times", colour="black",size=15),
        legend.title = element_blank())+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.05))+
  scale_x_continuous(expand = c(0, 0),limits = c(0,8.5),breaks = c(2,4,6,8))

ggsave("/mnt/disk6/3/3w_manuscript/fig2/distance/HiPore-C_all_dis.pdf", plot = p2, dpi = 300,width = 8,height = 6.5)
rm(H1,H2,H3,H4)

a1 <- fread("/mnt/disk1/1/gy/Cwalks/final_new_shuf/shuf_2w/hg38_cis_Cwalks_N4.txt2",header = F)
a2 <- fread("/mnt/disk1/1/gy/Cwalks/final_new_shuf/shuf_2w/hg38_cis_Cwalks_N4_3.txt2",header = F)
head(N_N_1)
N_N_3 <- a1 %>% inner_join(a2,by=c("V1"="V1","V4"="V4"))
N_N_3$dis <- abs(N_N_3$V2.x - N_N_3$V2.y)

#N_N_1 <- N_N_1 %>% sample_n(100000)
#N_N_2 <- N_N_2 %>% sample_n(100000)
#N_N_3 <- N_N_3 %>% sample_n(100000)
#N_N_k <- N_N_k %>% sample_n(100000)

merged <- c(N_N_1$dis,N_N_2$dis,N_N_3$dis,N_N_k$dis)
merged_final <- data.frame(MergedColumn = merged, DataFrame = rep(c("N/N+1","N/N+2","N/N+3","N/N+k"),
                                                                  times = c(nrow(N_N_1),nrow(N_N_2),nrow(N_N_3),nrow(N_N_k))))
head(merged_final)
colnames(merged_final) <- c("dis","type")
p2 <- ggplot(merged_final,aes(x = log10(dis),color=type))+
  stat_ecdf(geom="smooth",se=F,linewidth=3)+
  scale_color_manual(values = c('#386CAF','#7FC87F','#fc9533','#9f79d3'))+
  ylab("Cumulative frequency")+
  xlab("log10(distance(bp))")+
  theme_bw()+
  geom_hline(yintercept = 0.5,linetype = "dashed", color = "black")+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.ticks.x = element_line(colour = "black", size = 1),
        #axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.text.y=element_text(family="Times",size=17,colour="black"),
        axis.text.x=element_text(family="Times",size=17,colour="black"),
        axis.title.x=element_text(family="Times",size = 17,colour="black"),
        axis.title.y=element_text(family="Times",size = 17,colour="black"),
        legend.text=element_text(family="Times", colour="black",size=15),
        legend.title = element_blank())+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.05))+
  scale_x_continuous(expand = c(0, 0),limits = c(0,8.5),breaks = c(2,4,6,8))

ggsave("/mnt/disk6/3/3w_manuscript/fig2/distance/C-walks_all_dis.pdf", plot = p2, dpi = 300,width = 8,height = 6.5)




















