##=================================figS5b========================================================
a1 <- fread("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/chr16_bin_H3K27ac.tab",header = F)
a2 <- fread("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/chr16_bin_H3K27me3.tab",header = F)
a3 <- fread("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/chr16_bin_H3K9me3.tab",header = F)
a4 <- fread("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/chr16_bin_H3K4me3.tab",header = F)
a5 <- fread("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/chr16_bin_H3K36me3.tab",header = F)

A1 <- A1 %>% left_join(a1[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a2[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a3[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a4[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a5[,c(1,5)],by=c("V1"="V1")) %>% left_join(dd[,c(1,3)],by=c("V1"="V1"))

colnames(A1)[5:9] <- c("H3K27ac","H3K27me3","H3K9me3","H3K4me3","H3K36me3")

B1 <- B1 %>% left_join(a1[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a2[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a3[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a4[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a5[,c(1,5)],by=c("V1"="V1")) %>% left_join(dd[,c(1,3)],by=c("V1"="V1"))
colnames(B1)[5:9] <- c("H3K27ac","H3K27me3","H3K9me3","H3K4me3","H3K36me3")

B2 <- B2 %>% left_join(a1[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a2[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a3[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a4[,c(1,5)],by=c("V1"="V1")) %>% 
  left_join(a5[,c(1,5)],by=c("V1"="V1")) %>% left_join(dd[,c(1,3)],by=c("V1"="V1"))
colnames(B2)[5:9] <- c("H3K27ac","H3K27me3","H3K9me3","H3K4me3","H3K36me3")

Total <- rbind(A1,B1,B2)
Total$type <- c(rep("A1",nrow(A1)),rep("B1",nrow(B1)),rep("B2",nrow(B2)))

stats_1 <- Total %>% group_by(type) %>% summarise(median=median(H3K27ac),mean=mean(H3K27ac))
stats_2 <- Total %>% group_by(type) %>% summarise(median=median(H3K27me3),mean=mean(H3K27me3))
stats_3 <- Total %>% group_by(type) %>% summarise(median=median(H3K9me3),mean=mean(H3K9me3))
stats_4 <- Total %>% group_by(type) %>% summarise(median=median(H3K4me3),mean=mean(H3K4me3))
stats_5 <- Total %>% group_by(type) %>% summarise(median=median(H3K36me3),mean=mean(H3K36me3))


Total_stats <- rbind(stats_1,stats_2,stats_3,stats_4,stats_5)
Total_stats$chip <- c(rep("H3K27ac",3),rep("H3K27me3",3),rep("H3K9me3",3),rep("H3K4me3",3),rep("H3K36me3",3))

Total_stats$mean_chip <- c(rep(mean(a1$V5),3),rep(mean(a2$V5),3),rep(mean(a3$V5),3),rep(mean(a4$V5),3),rep(mean(a5$V5),3))

Total_stats$lfd <- Total_stats$mean/Total_stats$mean_chip

Total_stats_1 <- Total_stats[,c(1,6,4)]

heatmap_data <- reshape2::acast(Total_stats_1, type ~ chip, value.var = "lfd")
matrix_data <- as.matrix(heatmap_data)

bk <- c(seq(0,1,by=0.01),seq(1.01,2,by=0.01))

pdf('/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/ChIP_fd_rm0.2.pdf', width=5.8, height=3.5)
pheatmap(
  mat = matrix_data,
  color = c(colorRampPalette(colors = c("#2166AC","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B2182B"))(length(bk)/2)),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.1f",
  number_color = "black",
  fontsize_row = 15,
  fontsize_col = 15,
  border_color = "grey",
  cellwidth = 40,
  cellheight = 40,
  breaks = bk,fontsize=14
)
dev.off()
##========input AB infomation=========
#/mnt/disk5/2/3w3e/new_total/comp_bin/250Kb/PC1_AB_noNA.txt
#/mnt/disk5/2/3w3e/new_total/comp_bin/500Kb/PC1_AB_noNA.txt
ABcomp <- read.table("/mnt/disk5/2/3w3e/new_total/comp_bin/250Kb/PC1_AB_noNA.txt",header = F)

ABcomp$binID <- as.numeric(substring(ABcomp$V5,3))
nrow(ABcomp[ABcomp$V1=="chr16",])

PC1_matrix <- read.table("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/chr16_1D.matrix.txt",header = F)
PC1_matrix <- as.matrix(PC1_matrix)

library("factoextra")
fviz_nbclust(PC1_matrix, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2)

set.seed(123)
km_result <- kmeans(PC1_matrix,3,nstart = 24)
aa <- km_result$centers

pc1_pca_values <- fread("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/chr16/chr16_1D_PC1.txt",header = F)
dd <- cbind(pc1_pca_values,cluster = km_result$cluster)
dd <- dd %>% left_join(ABcomp[,c(1:3,5)],by=c("V1"="V5"))

#write.table(dd,file = "/mnt/disk6/3/3w_PCA/2D_PCA/chr18_kmeans_cluter4",col.names = T,row.names = F,quote = F,sep = "\t")
p1 <- fviz_cluster(km_result, data = PC1_matrix,
             palette = c("#7FC87F","#fc9533","#9f79d3"), ##chr16
             #palette = c("#9f79d3","#fc9533","#7FC87F"), ##chr19
             #palette = c("#fc9533","#7FC87F","#9f79d3"),
             #ellipse.level=7,
             ellipse.type = "none",
             star.plot = F,
             repel = F,
             ggtheme = theme_bw(),
             #ellipse.alpha=0.9,
             geom=c("point"),pointsize=0.8
)+scale_shape_manual(values=c(19,19,19)) ##c(15,19,1)
p1
 ggsave("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/Kmeans_plot_AB1B2.pdf", plot = p1, dpi = 600,width = 3.6,height = 3)

##different cluster  /mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr19/Kmeans
A1 <- dd[dd$cluster==1,c(4:6,1)]
B1 <- dd[dd$cluster==2,c(4:6,1)]
B2 <- dd[dd$cluster==3,c(4:6,1)]


##========================================fig5b=================================================
matrix_500Kb <- read.table("/mnt/disk6/3/3w_manuscript/all_m3/250kb_matrix3/OE/chr16_OE_250kb.matrix3",header = F,sep = "\t")
#matrix_500Kb_raw <- fread("/mnt/disk6/3/3w_manuscript/all_m3/250kb_matrix3/raw_count/chr19_250kb.matrix3",header = F)
#matrix_500Kb <- matrix_500Kb %>% left_join(matrix_500Kb_raw,by=c("V1"="V1","V2"="V2","V3"="V3"))
#matrix_500Kb <- matrix_500Kb[matrix_500Kb$V4.y>2,1:4]
#colnames(matrix_500Kb) <- c("V1","V2","V3","V4")
matrix_500Kb_1 <- matrix_500Kb[,c(1,3,2,4)]
matrix_500Kb_2 <- matrix_500Kb[,c(2,1,3,4)]
matrix_500Kb_3 <- matrix_500Kb[,c(2,3,1,4)]
matrix_500Kb_4 <- matrix_500Kb[,c(3,1,2,4)]
matrix_500Kb_5 <- matrix_500Kb[,c(3,2,1,4)]
colnames(matrix_500Kb_1) <- colnames(matrix_500Kb)
colnames(matrix_500Kb_2) <- colnames(matrix_500Kb)
colnames(matrix_500Kb_3) <- colnames(matrix_500Kb)
colnames(matrix_500Kb_4) <- colnames(matrix_500Kb)
colnames(matrix_500Kb_5) <- colnames(matrix_500Kb)

total_matrix_500Kb <- rbind(matrix_500Kb,matrix_500Kb_1,matrix_500Kb_2,matrix_500Kb_3,matrix_500Kb_4,matrix_500Kb_5)
rm(matrix_500Kb_1,matrix_500Kb_2,matrix_500Kb_3,matrix_500Kb_4,matrix_500Kb_5)

bed250Kb <- fread("/mnt/disk6/3/3w_manuscript/all_m3/250kb_matrix3/250000_abs.bed",header = F)
bed250Kb <- bed250Kb %>% inner_join(dd[,3:6],by=c("V1"="V1.y","V2"="V2.y","V3"="V3"))

matrix_500Kb_1 <- total_matrix_500Kb %>% inner_join(bed250Kb[,c(4,5)],by=c("V1"="V4")) %>% inner_join(bed250Kb[,c(4,5)],by=c("V2"="V4")) %>%
  inner_join(bed250Kb[,c(4,5)],by=c("V3"="V4"))

stats_matrix_500Kb_1 <- matrix_500Kb_1 %>% group_by(V3,cluster.x,cluster.y) %>% summarise(mean=mean(V4))
stats_matrix_500Kb_1 <- stats_matrix_500Kb_1 %>% left_join(bed250Kb[,c(4,5)],by=c("V3"="V4"))
stats_matrix_500Kb_1 <- stats_matrix_500Kb_1[stats_matrix_500Kb_1$cluster.x==stats_matrix_500Kb_1$cluster.y,]

A1 <- stats_matrix_500Kb_1[stats_matrix_500Kb_1$cluster.x==1,]
B1 <- stats_matrix_500Kb_1[stats_matrix_500Kb_1$cluster.x==2,]
B2 <- stats_matrix_500Kb_1[stats_matrix_500Kb_1$cluster.x==3,]

A1$mean <- (A1$mean - mean(A1$mean))/sd(A1$mean)
B1$mean <- (B1$mean - mean(B1$mean))/sd(B1$mean)
B2$mean <- (B2$mean - mean(B2$mean))/sd(B2$mean)

PIC_data <-rbind(A1[,c(1,4,5)],B1[,c(1,4,5)],B2[,c(1,4,5)])
PIC_data$d2 <- c(rep("A1_A1",nrow(A1)),rep("B1_B1",nrow(B1)),rep("B2_B2",nrow(B2)))

aaaaa <- unique(sort(PIC_data$d2))
heatmap_data_all <- matrix(0, nrow = 3, ncol = nrow(bed250Kb),dimnames = list(aaaaa,bed250Kb$V4))

rows <- match(PIC_data$d2,aaaaa)
cols <- match(PIC_data$V3,bed250Kb$V4)

heatmap_data_all[cbind(rows,cols)] <- PIC_data$mean


col_annot <- data.frame(Cluster = PIC_data$cluster[match(colnames(heatmap_data_all), PIC_data$V3)])
rownames(col_annot) <- colnames(heatmap_data_all)
col_annot$Cluster <- as.factor(col_annot$Cluster)
cluster_colors <- list(Cluster = c("1"="#7FC87F","2"="#fc9533","3" ="#9f79d3"))  ##chr16
#cluster_colors <- list(Cluster = c("1"="#9f79d3","2"="#fc9533","3" ="#7FC87F"))  ##chr19 chr15
#cluster_colors <- list(Cluster = c("1"="#fc9533","2"="#7FC87F","3" ="#9f79d3"))  ##chr17

#bk <- c(seq(0.3,2,by=0.01))
pdf('/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/2D_APA_cluster_filter0.2_scale.pdf', width=10, height=3)
pheatmap(
  mat = heatmap_data_all,
  color = colorRampPalette(colors = c("#2166AC","white","#B2182B"))(100),
  breaks = seq(from=-4, to=4, length.out=101),
  annotation_col = col_annot,
  annotation_colors = cluster_colors,
  cluster_rows = F,
  cluster_cols = F,
  scale = "none",
  show_rownames = T,
  fontsize_col = 8,
  border_color = NA,show_colnames=F
)

dev.off()


#/mnt/disk6/3/useful_file/K562_microC/K562_MicroC_R1R2R3R4_OE_250kb.txt
#/mnt/disk6/3/3w_manuscript/all_m3/split_3_2w/Total.split_3_2w_chr16_OE.txt
#/mnt/disk6/3/useful_file/K562_Hi-C/K562_MboI_total_chr16_OE_250kb.txt
w2_data <- fread("/mnt/disk6/3/useful_file/K562_Hi-C/K562_MboI_total_chr16_OE_250kb.txt",header = F)
w2_data_1 <- w2_data[,c(2,1,3)]
colnames(w2_data_1) <- colnames(w2_data)
Total_w2_data <- rbind(w2_data,w2_data_1)
rm(w2_data_1)

#chr16 10236 chr19 11512 chr17 10598 chr15 9828
Total_w2_data$V1 <- Total_w2_data$V1/250000 + 10236
Total_w2_data$V2 <- Total_w2_data$V2/250000 + 10236

w2_data_1 <- Total_w2_data %>% inner_join(bed250Kb[,4:5],by=c("V1"="V4"))
stats_w2_data_1 <- w2_data_1 %>% group_by(cluster,V2) %>% summarise(mean=mean(V3))

stats_w2_data_1 <- stats_w2_data_1[stats_w2_data_1$V2 %in% bed250Kb$V4,]
stats_w2_data_1 <- stats_w2_data_1 %>% group_by(cluster) %>% mutate(Zmeans=(mean-mean(mean))/sd(mean))
aaaaa <- unique(sort(stats_w2_data_1$cluster))
heatmap_data_2w_all <- matrix(0, nrow = 3, ncol = nrow(bed250Kb),dimnames = list(aaaaa,bed250Kb$V4))

rows <- match(stats_w2_data_1$cluster,aaaaa)
cols <- match(stats_w2_data_1$V2,bed250Kb$V4)


heatmap_data_2w_all[cbind(rows,cols)] <- stats_w2_data_1$Zmeans


pdf('/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/Hi-C_3w_split_3_APA_cluster_scale.pdf', width=10, height=3)
pheatmap(
  mat = heatmap_data_2w_all,
  color = colorRampPalette(colors = c("#2166AC","white","#B2182B"))(100),#color = colorRampPalette(colors = c("#2166AC","white","#B2182B"))(100),
  breaks = seq(from=-4, to=4, length.out=100),
  annotation_col = col_annot,
  annotation_colors = cluster_colors,
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  fontsize_col = 8,
  scale = "none",
  border_color = NA,show_colnames=F,
)

dev.off()



##======================fig5c===========================================
A1A1 <- PIC_data[PIC_data$d2=="B1_B1",]
W2_A1 <- stats_w2_data_1[stats_w2_data_1$cluster==2,]

A1A1 <- A1A1 %>% left_join(W2_A1[,c(2,4)],by=c("V3"="V2"))
cor.test(A1A1$mean,A1A1$Zmeans,method = "spearman")

A1A1$type <- ifelse(A1A1$mean*A1A1$Zmeans<0 & A1A1$mean>0,"diff1",
                    ifelse(A1A1$mean*A1A1$Zmeans<0 & A1A1$mean<0,"diff2",
                           ifelse(A1A1$mean*A1A1$Zmeans>0 & A1A1$mean<0,"left_down","right_up")))


stats_A1A1 <- A1A1 %>% group_by(cluster,type) %>% summarise(n=n())
stats_A1A1 <- stats_A1A1 %>% group_by(cluster) %>% mutate(ratio=n/sum(n))


A1A1$cluster <- as.factor(A1A1$cluster)
p1 <- ggplot(A1A1, aes(mean,Zmeans))+
  theme_classic()+
  geom_point(aes(color = cluster),size=0.5)+
  scale_color_manual(values = c("#7FC87F","#fc9533","#9f79d3")) +
  geom_point(
    data = subset(A1A1, cluster == "2"),aes(x = mean, y = Zmeans),color = "#fc9533",size = 0.5,alpha = 1) +
  xlim(-1.5,2)+ylim(-1.5,2)+    #B1 (-1.3,2)  #B2 -1.8,3 #chr17 -2,3
  labs(title = NULL, x = "B1B1 3w OE", y = "2w OE") + #distal proximal
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.line = element_line(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, colour = 'black'),
        axis.text.y = element_text(size = 15, colour = 'black'),
        #legend.position = c(0.7, 0.3),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, colour = 'black'))+
  #annotate("text", x = 1.5, y = 0.2, label = paste(round(nrow(cis[cis$E1>0 & cis$out5>0,]) / nrow(cis) * 100,1),"%",sep = ""), size = 5, colour = "red")+
  #annotate("text", x = -1.5, y = -0.2, label = paste(round(nrow(cis[cis$E1<0 & cis$out5<0,]) / nrow(cis) * 100,1),"%",sep = ""), size = 5, colour = "red")+
  #annotate("text", x = 2, y = 1, label = "R = 0.84", size = 5, colour = "black")+
  #geom_segment(aes(x = -2, y = -2, xend = 2, yend = 2), colour = "grey", size = 0.3)
geom_segment(aes(x = -1.5, y = 0, xend = 2, yend = 0), colour = "grey", size = 0.3)+
geom_segment(aes(x = 0, y = -1.5, xend = 0, yend = 2), colour = "grey", size = 0.3)
p1
ggsave("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/cor_OE_B1B1_insitu.pdf", plot = p1, dpi = 600,width = 4.7,height = 4.1)

##======================figS5c 3w OE 3w_2D OE and in situ OE=============================

w2_data <- fread("/mnt/disk6/3/useful_file/K562_Hi-C/K562_MboI_total_chr16_OE_250kb.txt",header = F)
w2_data_2 <- w2_data[,c(2,1,3)]
colnames(w2_data_2) <- colnames(w2_data)
Total_w2_data <- rbind(w2_data,w2_data_2)
rm(w2_data_2)

#chr16 10236 chr19 11512 chr17 10598 chr15 9828
Total_w2_data$V1 <- Total_w2_data$V1/250000 + 10236
Total_w2_data$V2 <- Total_w2_data$V2/250000 + 10236

w2_data_1 <- Total_w2_data %>% inner_join(bed250Kb[,4:5],by=c("V1"="V4"))
stats_w2_data_2 <- w2_data_1 %>% group_by(cluster,V2) %>% summarise(mean=mean(V3))

stats_w2_data_2 <- stats_w2_data_2[stats_w2_data_2$V2 %in% bed250Kb$V4,]
stats_w2_data_2 <- stats_w2_data_2 %>% group_by(cluster) %>% mutate(Zmeans=(mean-mean(mean))/sd(mean))

A1A1 <- PIC_data[PIC_data$d2=="A1_A1",]
W2_A1 <- stats_w2_data_1[stats_w2_data_1$cluster==1,]
IN_A1 <- stats_w2_data_2[stats_w2_data_2$cluster==1,]

Total_matrix <- A1A1 %>% left_join(W2_A1[,c(2,4)],by=c("V3"="V2")) %>% left_join(IN_A1[,c(2,4)],by=c("V3"="V2"))

Total_matrix <- Total_matrix[,c(2,5,6)]
colnames(Total_matrix) <- c("3wHi-C","3wHi-C 2D","in situ Hi-C")

corr <- round(cor(Total_matrix),2)

#color_palette <- colorRampPalette(c("blue","red"))(50)
pdf("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/A1A1_3_method_heatmap.pdf",width=4.6, height=3.8, bg = "transparent")
pheatmap(corr, cluster_cols=F, cluster_rows=F,
         legend = T,
         border_color=NA,
         show_rownames = T, show_colnames = T,
         color = colorRampPalette(colors = c("white","#FF6347"))(100),
         breaks = seq(0.65,1,length.out = 100),
         legend_breaks = seq(0.6,1,0.1),
         fontsize = 18,
         display_numbers = TRUE,
         number_format = "%.2f")
dev.off()


##======================fig5d===========================================
#/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans
segregated_data <- fread("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/B1B2_segregated.txt",header=F)

segregated_data_1 <- segregated_data[segregated_data$V5=="B2-B2",]

blank_data <- data.frame(group = c("B2","B2","A1","A1","A1","A1","A1"),x = 0, y = c(0,0.8,0,0.2,0.4,0.81,0.89))

p1 <- ggplot() +
  geom_col(data = segregated_data_1,position = position_dodge(width = 0.8), width = 0.7, aes(x = V3, y = V1, fill = V2)) +
  #geom_blank(data = blank_data, aes(x = x, y = y)) +
  facet_wrap(~ V4, scales = "free_y", ncol = 1) +
  scale_y_continuous(limits = c(0,0.9),breaks = c(0,0.2,0.4,0.6,0.8))+
  scale_fill_manual(values = c('#386CAF','#9f79d3','#fc9533','grey','blue','red','black','#7FC87F'))+
  labs(y = "Percentage",x="") +
  theme_bw()+
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
ggsave("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/B2B2-B1A1_segregated_all.pdf", plot = p1, dpi = 300,width = 4.3,height = 4.3)


##==========================================fig5e and figS5e======================================

dd$binID <- substring(dd$V1,3)
dd$binID <- as.numeric(dd$binID)
head(B1B1B2)
B1 <- c(dd[dd$cluster==2,]$binID)
B2 <- c(dd[dd$cluster==3,]$binID)

B1B1B2 <- matrix_500Kb_1[matrix_500Kb_1$cluster.x==3 & matrix_500Kb_1$cluster.y==3,]

##2w
w2_data_1_1 <- w2_data_1 %>% inner_join(bed250Kb[,4:5],by=c("V2"="V4"))

for (i in 33:40) {
  ##3w
  for (j in i+1:41) {
    B1B1B2_1 <- B1B1B2[B1B1B2$V1==B2[i] & B1B1B2$V2==B2[j],]
    B1B1B2_1$V4 <- (B1B1B2_1$V4 - mean(B1B1B2_1$V4))/sd(B1B1B2_1$V4)
    B1B1B2_1$type <- "B2-B2"
    ##2w
    w2_data_1_1_1 <- w2_data_1_1[w2_data_1_1$V1==B2[i],]
    w2_data_1_1_1$V3 <- (w2_data_1_1_1$V3 - mean(w2_data_1_1_1$V3))/sd(w2_data_1_1_1$V3)
    
    overlap_1 <- data.frame(bin=intersect(B1B1B2_1$V3,w2_data_1_1_1$V2))
    overlap_1 <- overlap_1 %>% left_join(bed250Kb[,4:5],by=c("bin"="V4"))
    overlap_1 <- overlap_1[order(overlap_1$bin),]
    aaaaa <- "B2-B2"
    heatmap_data_3w_all <- matrix(0, nrow = 1, ncol = nrow(overlap_1),dimnames = list(aaaaa,overlap_1$bin))
    
    rows <- match(B1B1B2_1$type,aaaaa)
    cols <- match(B1B1B2_1$V3,overlap_1$bin)
    
    heatmap_data_3w_all[cbind(rows,cols)] <- B1B1B2_1$V4
    
    head(overlap_1)
    
    col_annot <- data.frame(Cluster = overlap_1$cluster[match(colnames(heatmap_data_3w_all),overlap_1$bin)])
    rownames(col_annot) <- colnames(heatmap_data_3w_all)
    col_annot$Cluster <- as.factor(col_annot$Cluster)
    cluster_colors <- list(Cluster = c("1"="#7FC87F","2"="#fc9533","3" ="#9f79d3"))
    
    w2_data_1_1_1$type <- "B2"
    bbbb <- "B2"
    w2_data_1_1_1 <- w2_data_1_1_1[w2_data_1_1_1$V2 %in% overlap_1$bin,]
    heatmap_data_2w_all <- matrix(0, nrow = 1, ncol = nrow(overlap_1),dimnames = list(bbbb,overlap_1$bin))
    
    rows <- match(w2_data_1_1_1$type,bbbb)
    cols <- match(w2_data_1_1_1$V2,overlap_1$bin)
    
    heatmap_data_2w_all[cbind(rows,cols)] <- w2_data_1_1_1$V3
    
    Total_heatmap_2w_3w <- rbind(heatmap_data_3w_all,heatmap_data_2w_all)
    output_file <- paste("/mnt/disk2/3/gy/fig5_new/3w_穿透/250Kb/manuscript_related/chr16/Kmeans/all_B2B2_3W_2D/",B2[i],"_",B2[j],".pdf",sep = "")
    pdf(output_file, width=10, height=3)
    pheatmap(
      mat = Total_heatmap_2w_3w,
      color = colorRampPalette(colors = c("#2166AC","white","#B2182B"))(100),#color = colorRampPalette(colors = c("#2166AC","white","#B2182B"))(100),
      breaks = seq(from=-4, to=4, length.out=100),
      annotation_col = col_annot,
      annotation_colors = cluster_colors,
      cluster_rows = F,
      cluster_cols = F,
      show_rownames = T,
      fontsize_col = 8,
      scale = "none",
      border_color = NA,show_colnames=F,
    )
    dev.off()
  }
}


