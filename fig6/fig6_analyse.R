##===========================================fig6c=============================================

##using END-seq peak: APH specific
DSB_peak <- fread("APH_s_10Kb_bin.txt",header = F)
DSB_peak <- DSB_peak[DSB_peak$V1=="chr1",]

##input 3wHi-C data
APH_3w <- fread("/mnt/disk2/3/gy/APH_DMSO/APH/10kb_m3/split/APH_chr1_intra_w3_10kb.matrix3",header = F)
DMSO_3w <- fread("/mnt/disk2/3/gy/APH_DMSO/DMSO/10kb_m3/split/DMSO_chr1_intra_w3_10kb.matrix3",header = F)

##count normalized
APH_3w$V4 <- APH_3w$V4/sum(APH_3w$V4)*1000000
DMSO_3w$V4 <- DMSO_3w$V4/sum(DMSO_3w$V4)*1000000

Total_3w <- APH_3w %>% inner_join(DMSO_3w,by=c("V1"="V1","V2"="V2","V3"="V3"))
colnames(Total_3w)[4:5] <- c("APH","DMSO")

Total_3w$lfd <- log2(Total_3w$APH/Total_3w$DMSO)

##all or within One TAD
bed10Kb <- fread("/mnt/disk2/3/gy/fig6/TAD_related_DSB/10Kb/10Kb_bin_overlap_TAD.txt",header = F)

Total_3w_1 <- Total_3w %>% inner_join(bed10Kb[,c(4,5)],by=c("V1"="V4")) %>% inner_join(bed10Kb[,c(4,5)],by=c("V2"="V4")) %>% inner_join(bed10Kb[,c(4,5)],by=c("V3"="V4"))
#Total_3w_1 <- Total_3w_1[Total_3w_1$V5.x==Total_3w_1$V5.y & Total_3w_1$V5.x==Total_3w_1$V5 & Total_3w_1$V5.y==Total_3w_1$V5,]

##DSB number class
Total_3w_3_DSB <- Total_3w_1[Total_3w_1$V1 %in% DSB_peak$V4 & Total_3w_1$V2 %in% DSB_peak$V4 & Total_3w_1$V3 %in% DSB_peak$V4,]
Total_3w_0_DSB <- Total_3w_1[!(Total_3w_1$V1 %in% DSB_peak$V4) & !(Total_3w_1$V2 %in% DSB_peak$V4) & !(Total_3w_1$V3 %in% DSB_peak$V4),]
Total_3w_2_DSB <- Total_3w_1[
  rowSums(
    cbind(
      Total_3w_1$V1 %in% DSB_peak$V4,
      Total_3w_1$V2 %in% DSB_peak$V4,
      Total_3w_1$V3 %in% DSB_peak$V4
    )
  ) == 2,
]

Total_3w_1_DSB <- Total_3w_1[
  rowSums(
    cbind(
      Total_3w_1$V1 %in% DSB_peak$V4,
      Total_3w_1$V2 %in% DSB_peak$V4,
      Total_3w_1$V3 %in% DSB_peak$V4
    )
  ) == 1,
]

##except DMSO specific sites
DMSO_S <- fread("DMSO_s_10Kb_bin.txt",header = F)
Total_3w_0_DSB <- Total_3w_0_DSB[!(Total_3w_0_DSB$V1 %in% DMSO_S$V4) & !(Total_3w_0_DSB$V2 %in% DMSO_S$V4) & !(Total_3w_0_DSB$V3 %in% DMSO_S$V4),]
Total_3w_1_DSB <- Total_3w_1_DSB[!(Total_3w_1_DSB$V1 %in% DMSO_S$V4) & !(Total_3w_1_DSB$V2 %in% DMSO_S$V4) & !(Total_3w_1_DSB$V3 %in% DMSO_S$V4),]
Total_3w_2_DSB <- Total_3w_2_DSB[!(Total_3w_2_DSB$V1 %in% DMSO_S$V4) & !(Total_3w_2_DSB$V2 %in% DMSO_S$V4) & !(Total_3w_2_DSB$V3 %in% DMSO_S$V4),]

Total_3w_0_3_DSB <- rbind(Total_3w_0_DSB,Total_3w_3_DSB)
Total_3w_0_3_DSB$type <- c(rep("Non_DSB",nrow(Total_3w_0_DSB)),rep("DSB",nrow(Total_3w_3_DSB)))

Total_3w_0_3_DSB$contact_type <- ifelse(Total_3w_0_3_DSB$lfd>log2(1.5),"3w_UP",
                                        ifelse(Total_3w_0_3_DSB$lfd < -log2(1.5),"3w_Down","others")) 

Total_3w_0_3_DSB$span <- Total_3w_0_3_DSB$V3 - Total_3w_0_3_DSB$V2

stats_Total_3w_0_3_DSB <- Total_3w_0_3_DSB %>% group_by(type,contact_type) %>% summarise(n=n())
stats_Total_3w_0_3_DSB <- stats_Total_3w_0_3_DSB %>% group_by(type) %>% mutate(ratio=n/sum(n))

stats_Total_3w_0_3_DSB_1 <- stats_Total_3w_0_3_DSB[stats_Total_3w_0_3_DSB$contact_type != "others",]

p1 <- ggplot(stats_Total_3w_0_3_DSB_1, aes(x = type, y = ratio, fill = contact_type)) +
  geom_col(stat = 'identity',position = "stack",color=NA, width=0.55) + 
  scale_fill_manual(values = c("#5479BD","#CD575B"))+
  labs(y = "Value",x="") +
  theme_bw()+
  ylim(0,0.55)+
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.line = element_line(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, colour = 'black'),
        axis.text.y = element_text(size = 15, colour = 'black'),
        legend.position = "top",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, colour = 'black'),
        strip.background = element_rect(fill = "white",linetype = "dotted"),
        strip.text = element_text(color = "black",
                                  face = "bold",
                                  size = 14))
p1
ggsave("/mnt/disk2/3/gy/fig6/TAD_related_DSB/10Kb/contact_up_down_percentage.pdf", plot = p1, dpi = 600,width = 3.2,height = 3.5)

##===========================================fig6d=============================================
#span distance
Total_3w_3_DSB_dis <- Total_3w_3_DSB
Total_3w_3_DSB_dis$dis <- Total_3w_3_DSB_dis$V3 - Total_3w_3_DSB_dis$V2
Total_3w_3_DSB_dis$contact_type <- ifelse(Total_3w_3_DSB_dis$lfd>log2(1.5),"3w_UP",
                                        ifelse(Total_3w_3_DSB_dis$lfd < -log2(1.5),"3w_Down","others")) 

Total_3w_3_DSB_dis$dis_type <- ifelse(Total_3w_3_DSB_dis$dis<=5,"<=5",
                                      ifelse(Total_3w_3_DSB_dis$dis>5 & Total_3w_3_DSB_dis$dis<=10,"5-10",ifelse(Total_3w_3_DSB_dis$dis>10 & Total_3w_3_DSB_dis$dis<=15,"10-15",">15")))

stats_Total_3w_3_DSB_dis <- Total_3w_3_DSB_dis %>% group_by(dis_type,contact_type) %>% summarise(n=n())
stats_Total_3w_3_DSB_dis <- stats_Total_3w_3_DSB_dis %>% group_by(dis_type) %>% mutate(ratio=n/sum(n))

stats_Total_3w_3_DSB_dis <- stats_Total_3w_3_DSB_dis[stats_Total_3w_3_DSB_dis$contact_type != "others",]

AA <- factor(stats_Total_3w_3_DSB_dis$dis_type,levels = c("<=5","5-10","10-15",">15"))
p1 <- ggplot(stats_Total_3w_3_DSB_dis, aes(x = AA, y = ratio, fill = contact_type)) +
  geom_col(stat = 'identity',position = "stack",color=NA, width=0.55) +
  scale_fill_manual(values = c("#5479BD","#CD575B"))+
  labs(y = "Value",x="") +
  theme_bw()+
  ylim(0,0.6)+
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.line = element_line(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, colour = 'black'),
        axis.text.y = element_text(size = 15, colour = 'black'),
        legend.position = "top",
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, colour = 'black'),
        strip.background = element_rect(fill = "white",linetype = "dotted"),
        strip.text = element_text(color = "black",
                                  face = "bold",
                                  size = 14))
p1
ggsave("/mnt/disk2/3/gy/fig6/TAD_related_DSB/10Kb/DSB_group_contact_up_down_dis_bar.pdf", plot = p1, dpi = 600,width = 4,height = 3.5)


##===========================================fig6g, figS6f=============================================
Total_pic2 <- rbind(Total_3w_1_DSB,Total_3w_2_DSB,Total_3w_3_DSB)

Total_pic2$classss <- c(rep("1_DSB",nrow(Total_3w_1_DSB)),rep("2_DSB",nrow(Total_3w_2_DSB)),rep("3_DSB",nrow(Total_3w_3_DSB)))
stats_Total_pic2 <- Total_pic2 %>% group_by(classss) %>% summarise(mean=mean(lfd),median=median(lfd))

##bar plot
p <- ggplot(stats_Total_pic2,aes(x=classss,y=mean,fill=classss))+
  geom_bar(stat = 'identity',position = "stack",color=NA, width=0.6)+
  scale_fill_manual(values = c('#386CAF','#fc9533','#9f79d3','#7FC87F','grey','blue','red','black'))+
  scale_color_manual(values = c("black", "black", "black")) +
  theme_classic()+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.5),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.text.y=element_text(size=17,colour="black"),
        axis.text.x=element_text(size=17,colour="black"),
        axis.title.x=element_text(size = 17,colour="black"),
        axis.title.y=element_text(size = 17,colour="black"),
        legend.position = "none",legend.text = element_text(size = 16,colour="black"),
        legend.title = element_blank())+
  xlab("")+ylab("Mean log2 (Fold change)")   #ylab("Sample (N)")
p

p1 <- ggplot(Total_pic2,aes(x = lfd,color=classss))+
  stat_ecdf(geom="smooth",se=F,linewidth=1.2)+
  scale_color_manual(values = c('#386CAF','#fc9533','#9f79d3','#7FC87F'))+
  ylab("Cumulative frequency")+
  xlab("lfd")+
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
        legend.title = element_blank())
p1
ggsave("/mnt/disk2/3/gy/fig6/TAD_related_DSB/10Kb/1_3DSB_within_TAD_lfd_cumulative.pdf", plot = p1, dpi = 600,width = 5,height = 3.5)


##===========================================fig6f=============================================
APH_w3_1 <- fread("/mnt/disk2/3/gy/APH_DMSO/APH/filter_gt500_APH_f_m3_Total_0.bed",header = F)
APH_w3_2 <- fread("/mnt/disk2/3/gy/APH_DMSO/APH/filter_gt500_APH_f_m3_Total_1.bed",header = F)
APH_w3_3 <- fread("/mnt/disk2/3/gy/APH_DMSO/APH/filter_gt500_APH_f_m3_Total_2.bed",header = F)

Total_w3 <- cbind(APH_w3_1[,c(1,3,4)],APH_w3_2[,c(1,3)],APH_w3_3[,c(1,3)])
rm(APH_w3_1,APH_w3_2,APH_w3_3)

colnames(Total_w3) <- paste("V",1:7,sep = "")

DSB_1 <- fread("APH_DSB_1",header = F)
DSB_2 <- fread("APH_DSB_2",header = F)
DSB_3 <- fread("APH_DSB_3",header = F)

DSB_1_3w <- Total_w3[Total_w3$V3 %in% DSB_1$V1,]
DSB_2_3w <- Total_w3[Total_w3$V3 %in% DSB_2$V1,]
DSB_3_3w <- Total_w3[Total_w3$V3 %in% DSB_3$V1,]

##calculate span distance
Total_DSB_3w_APH <- rbind(DSB_1_3w,DSB_2_3w,DSB_3_3w)
Total_DSB_3w_APH$type <- c(rep("DSB_1",nrow(DSB_1_3w)),rep("DSB_2",nrow(DSB_2_3w)),rep("DSB_3",nrow(DSB_3_3w)))

##cis
Total_DSB_3w_APH <- Total_DSB_3w_APH[Total_DSB_3w_APH$V1==Total_DSB_3w_APH$V4 & Total_DSB_3w_APH$V1==Total_DSB_3w_APH$V6 & Total_DSB_3w_APH$V4==Total_DSB_3w_APH$V6,]
Total_DSB_3w_APH <- Total_DSB_3w_APH[,c(2,5,7,8)]
Total_DSB_3w_APH[,1:3] <- data.frame(t(apply(Total_DSB_3w_APH[,1:3],1,function(x) sort(x))))
Total_DSB_3w_APH$span <- Total_DSB_3w_APH$V7 - Total_DSB_3w_APH$V2
total_APH_DMSO <- rbind(Total_DSB_3w_APH,Total_DSB_3w)
total_APH_DMSO$group <- c(rep("APH",nrow(Total_DSB_3w_APH)),rep("DMSO",nrow(Total_DSB_3w)))
p1 <- ggplot(total_APH_DMSO[total_APH_DMSO$group=="APH",],aes(x = log10(span),color=type))+
  stat_ecdf(geom="smooth",se=F,linewidth=1.1)+
  scale_color_manual(values = c('#386CAF','#fc9533','#9f79d3','#7FC87F'))+
  ylab("Cumulative frequency")+
  xlab("log10(distance(bp))")+
  theme_bw()+
  geom_hline(yintercept = 0.5,linetype = "dashed", color = "black")+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        #axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.text.y=element_text(family="Times",size=17,colour="black"),
        axis.text.x=element_text(family="Times",size=17,colour="black"),
        axis.title.x=element_text(family="Times",size = 17,colour="black"),
        axis.title.y=element_text(family="Times",size = 17,colour="black"),
        legend.text=element_text(family="Times", colour="black",size=15),
        legend.title = element_blank())
p1
ggsave("/mnt/disk2/3/gy/fig6/TAD_related_DSB/fragment_DSB/APH_1_3_DSB_3w_frag.tiff", plot = p1, dpi = 600,width = 5.1,height = 3.5)



