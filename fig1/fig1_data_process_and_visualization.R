##=================================fig1d: Multiple index vs. Non-multiple index pairing time====================================
##/mnt/disk6/3/3w_manuscript/fig1_script/d
data <- read.table("merge_CPUtime_final.txt",header = F)
data$V4 <- c(rep("multiplxing",18),rep("non-multiplxing",18))
stats <- data %>% group_by(V1,V4) %>% summarise(mean=mean(V3),sd=sd(V3))
stats$mean <- stats$mean/min(stats$mean)
pdf("multiplx_non_CPUtime_new.pdf", width = 4.5, height = 3.5, bg = "transparent")
ggplot(stats,aes(x=V1/1000000,y=mean,color=V4)) +
  geom_line(aes(color = V4),linewidth=0.8) +
  geom_point(size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, colour = "black", size = 15, angle = 0),
        axis.text.y = element_text(size = 16, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 16, face = "plain", colour = "black"),
        axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 2, hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.5),
        axis.line.x = element_line(colour = "black", size = 0.5),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "lines"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "plain", colour = "black"),
        plot.title = element_text(colour = "black", face = "plain", size = 15, vjust = 8, hjust = 1)) +
  scale_x_continuous(breaks = c(5,10,15,20,25,30))+
  labs(title = "", x = "Number of read pairs (M)", y = "Relative Pairing time")+
  scale_color_manual(values = c("#386CAF", "#DB9052"))
dev.off()

##===========================fig1e: The overall pairing efficiency of single-indexed 3wHi-C data with accumulation of reads===================================
##/mnt/disk6/3/3w_manuscript/质控/financy_curve
financy_curve <- read.table("financy.curve.txt",header = F)
colors <- colorRampPalette(c("#DB9052","#386CAF"))(96)
all_colors <- colors()
random_colors <- sample(all_colors, 96)

pdf("financy_curve2.pdf", width = 3.8, height = 3.5, bg = "transparent")
ggplot(data = financy_curve,mapping = aes(x=V2/1000000,y=V8,color=V1))+
  geom_line(linewidth=0.8)+
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, colour = "black", size = 15, angle = 0),
        axis.text.y = element_text(size = 16, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 16, face = "plain", colour = "black"),
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
        legend.position = "none",
        plot.title = element_text(colour = "black", face = "plain", size = 15, vjust = 8, hjust = 1)) +
  labs(title = "", x = "Raw reads (M)", y = "Efficiency\n(Contribution / Rawreads)")+
  scale_color_manual(values = random_colors)+
  scale_y_continuous(expand = c(0, 0),limits = c(0,0.058))+scale_x_continuous(limits = c(0,13),breaks=c(1,3,5,7,9,11))
dev.off()

##=================================fig1f: PCC of 3wHi-C rep1-4=============================================
##/mnt/disk6/3/3w_manuscript/all_m3/repli_4
R1 <- fread("R1_total_250kb.matrix3",header = F)
R2 <- fread("R2_total_250kb.matrix3",header = F)
R3 <- fread("R3_total_250kb.matrix3",header = F)
R4 <- fread("R4_total_250kb.matrix3",header = F)
inner_result <- R1 %>% inner_join(R2,by=c("V1"="V1","V2"="V2","V3"="V3")) %>%
  inner_join(R3,by=c("V1"="V1","V2"="V2","V3"="V3")) %>%
  inner_join(R4,by=c("V1"="V1","V2"="V2","V3"="V3"))

colnames(inner_result)[4:7] <- c("R1","R2","R3","R4")

corr <- round(cor(inner_result[,4:7]),2)
p.mat <- cor_pmat(inner_result[,4:7])
library("pheatmap")
pdf("r1_r4_3w_cor_2.pdf",width=5.1, height=4.2, bg = "transparent")
pheatmap(corr, cluster_cols=F, cluster_rows=F,
         legend = T,
         border_color=NA,
         show_rownames = T, show_colnames = T,
         color = colorRampPalette(colors = c("white","#FF6347"))(100),
         breaks = seq(0.6,1,length.out = 100),
         legend_breaks = seq(0.6,1,0.1),
         fontsize = 18,
         display_numbers = TRUE,
         number_format = "%.2f")
dev.off()

##=================================fig1i: Overlap Venn diagram=============================================
#The detailed code for the transitivity analysis can be found in the file "cliques.sh".
w3_Hi_C <- fread("/mnt/disk6/3/3w_manuscript/all_m3/100kb_matrix3/Total_intra_w3_raw_100000.matrix3",header = F)
w3_Hi_C <- w3_Hi_C[w3_Hi_C$V4>2,]

cliques_file <- read.table("/mnt/disk6/3/useful_file/scHi-C/GAGE-seq/total_cell_cliques_100kb_3w",header = F,sep = " ")
sorted_values <- apply(cliques_file[,2:4], 1, function(x) sort(x))
sorted_df <- data.frame(t(sorted_values))
cliques_file[,2:4] <- sorted_df
cliques_file[,2:4] <- cliques_file[,2:4] + 1
rm(sorted_values)

cell_2w <- fread("/mnt/disk6/3/useful_file/scHi-C/GAGE-seq/K562_GAGE_seq.freq_100kb.txt",header = F)

cliques_file_2w_f <- inner_join(cliques_file,cell_2w,by=c("V1"="V1","V2"="V2","V3"="V3"))
cliques_file_2w_f <- inner_join(cliques_file_2w_f,cell_2w,by=c("V1"="V1","V2"="V2","V4.x"="V3"))
cliques_file_2w_f <- inner_join(cliques_file_2w_f,cell_2w,by=c("V1"="V1","V3"="V2","V4.x"="V3"))

cliques_file_2w_f$min <- apply(cliques_file_2w_f[,5:7],1,function(x) min(x))
cliques_file_2w_final <- cliques_file_2w_f[cliques_file_2w_f$min>2,c(1:4,8)]

cell_num <- unique(sort(cliques_file_2w_final$V1))
cell_stats <- cliques_file_2w_final %>% group_by(V2,V3,V4.x) %>% summarise(count=dplyr::n(),sum=sum(min))

inner_result <- inner_join(cell_stats,w3_Hi_C,by=c("V1"="V1","V2"="V2","V3"="V3"))

`3wHi-C` <- 1:3270397
`GAGE-seq` <- 3236724:3278248

s1 <- list(`3wHi-C` = `3wHi-C`, 
           `GAGE-seq` = `GAGE-seq`)

v1 <- venn.diagram(x = s1,filename = NULL,
                   resolution = 500,
                   scaled = T,
                   alpha=c(0.6,0.6),
                   lwd = 1,
                   lty = 1,
                   col = 'black',
                   fill = c("#C36518","#a085bd"),
                   cex = 1.2,
                   fontface = "bold",
                   fonrfamily = "sans",
                   cat.cex = 1.3,
                   cat.col = "black",
                   cat.fontface = "bold",
                   cat.default.pos = "outer",
                   cat.pos = c(-10,20),
                   cat.dist = c(0.1,0.1),
                   cat.fontfamily = "sans")
pdf("Venn_compare_100kb_GAGEseq_gt2_cell.pdf", width = 3.2, height = 3.1, bg = "transparent")
cowplot::plot_grid(v1)
dev.off()

##===============================fig1j: The mean spatial distances of three-way contacts with different contact frequencies===============================
imaging <- fread("/mnt/disk6/3/software/ChromatinImaging-master/Data/probe_distances_parallel.csv2",header = F)
head(imaging)
imaging[,2:3] <- data.frame(t(apply(imaging[,2:3], 1, function(x) sort(x))))

imaging_stats <- imaging %>%
  group_by(V2,V3) %>%
  slice_min(V4,n=10) %>%
  summarise(mean_top100 = mean(V4))

my_3w <- read.table("/mnt/disk6/3/software/ChromatinImaging-master/Data/mydata/chr21_3w.result.txt",header = F)
head(my_3w)
my_3w[,1:3] <- data.frame(t(apply(my_3w[,1:3], 1, function(x) sort(x))))

stats_my_3w <- my_3w %>% group_by(V1,V2,V3) %>% summarise(count=dplyr::n())

stats_my_3w_1 <- stats_my_3w[stats_my_3w$V1 != stats_my_3w$V2 & stats_my_3w$V2 != stats_my_3w$V3 & stats_my_3w$V1 != stats_my_3w$V3,]
stats_my_3w_1[,1:3] <- stats_my_3w_1[,1:3] + 1

stats_my_3w_1 <- left_join(stats_my_3w_1,imaging_stats,by=c("V1"="V2","V2"="V3"))
stats_my_3w_1 <- left_join(stats_my_3w_1,imaging_stats,by=c("V1"="V2","V3"="V3"))
stats_my_3w_1 <- left_join(stats_my_3w_1,imaging_stats,by=c("V2"="V2","V3"="V3"))

stats_my_3w_1$mean_dis <- apply(stats_my_3w_1[,5:7],1,function(x) mean(x))

stats_my_3w_1$Type <- ifelse(stats_my_3w_1$count>=1 & stats_my_3w_1$count<=2,"[1,2]",
                             ifelse(stats_my_3w_1$count>2 & stats_my_3w_1$count<=4,"(2,4]",
                                    ifelse(stats_my_3w_1$count>4 & stats_my_3w_1$count<=6,"(4,6]","(6,)")))

aa <- factor(stats_my_3w_1$Type,levels = c("[1,2]","(2,4]","(4,6]","(6,)"))
pdf("/mnt/disk6/3/3w_manuscript/fig1_script/h/distance_VS_raw_count_imaging_top10.pdf", width = 4.2, height = 3.2, bg = "transparent")
ggplot(data = stats_my_3w_1,mapping = aes(x=aa,y=mean_dis))+
  geom_boxplot(outlier.colour = NA,fill="#B9B3E0")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.line = element_line(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, colour = 'black'),
        axis.text.y = element_text(size = 15, colour = 'black'),
        legend.position = "none")+xlab("Raw count")+ylab("Mean spatial distance (nm)") + coord_cartesian(ylim = c(0,45))
dev.off()











