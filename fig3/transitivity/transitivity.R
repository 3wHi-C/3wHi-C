setwd("/mnt/disk2/3/gy/3w_fig3/cliques")

##chr1
#bed_50kb <- fread("/mnt/disk6/3/useful_file/K562_microC/matrix/K562_MicroC_R1+R2+R3+R4.50000_abs.bed",header = F)
#bed_50kb_stats <- bed_50kb %>% group_by(V1) %>% summarise(min=min(V4),max=max(V4))

matrix_50kb <- fread("/mnt/disk6/3/useful_file/K562_microC/matrix/K562_MicroC_R1+R2+R3+R4.50000_iced.matrix",header = F)
matrix3_50kb <- fread("/mnt/disk6/3/3w_manuscript/all_m3/50kb_matrix3/Total_w3_intra_50000.matrix3",header = F)
matrix3_50kb <- matrix3_50kb[matrix3_50kb$V1<=4980 & matrix3_50kb$V2<=4980 & matrix3_50kb$V3<=4980 & matrix3_50kb$V4>2,]

##filter chr1
matrix_50kb <- matrix_50kb[matrix_50kb$V1<=4980 & matrix_50kb$V2<=4980,]

##cliques
matrix_50kb <- matrix_50kb[matrix_50kb$V1 != matrix_50kb$V2 & matrix_50kb$V3>quantile(matrix_50kb$V3,0.9),]

##step 1
matrix_50kb_1 <- inner_join(matrix_50kb,matrix_50kb,by=c("V1"="V1"),relationship = "many-to-many")
matrix_50kb_1 <- matrix_50kb_1[matrix_50kb_1$V2.x != matrix_50kb_1$V2.y,]

##step 2
matrix_50kb_1 <- inner_join(matrix_50kb_1,matrix_50kb,by=c("V2.x"="V1","V2.y"="V2"))
matrix_50kb_1 <- matrix_50kb_1[,c(1,2,4,3,5,6)]
colnames(matrix_50kb_1) <- c("V1","V2","V3","Count12","Count13","Count23")

matrix_50kb_1$min <- apply(matrix_50kb_1[,4:6],1,function(x) min(x))

##互作距离分类
##GAM: We then selected the strongest 2% of all Hi-C triplets from every genomic distance 
##for which there were at least 500 possible triplets
matrix_50kb_1$dis1 <- abs(matrix_50kb_1$V2 - matrix_50kb_1$V1)
matrix_50kb_1$dis2 <- abs(matrix_50kb_1$V3 - matrix_50kb_1$V2)

dis_stats <- matrix_50kb_1 %>% group_by(dis1,dis2) %>% summarise(n=n())
dis_stats_1 <- dis_stats[dis_stats$n>=500,]
head(matrix_50kb_2)

matrix_50kb_2 <- inner_join(matrix_50kb_1,dis_stats_1,by=c("dis1"="dis1","dis2"="dis2"))

##根据每个互作距离挑选前2%的互作位点作为三元互作
head(total_compare)
total_compare <- matrix_50kb_2 %>% group_by(dis1,dis2) %>%
  slice_max(min,prop = 0.02)

##计算两者的overlap
overlap_data <- inner_join(total_compare[,c(1,2,3,7)],matrix3_50kb,by=c("V1"="V1","V2"="V2","V3"="V3"))

##绘制Veen图
`3wHi-C` <- 1:358559
`Predicted 3-ways` <- 344710:684944
s1 <- list(`3wHi-C` = `3wHi-C`, 
           `Predicted` = `Predicted 3-ways`)
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
pdf("/mnt/disk2/3/gy/3w_fig3/cliques/Venn_compare_50kb.pdf", width = 3.2, height = 3.1, bg = "transparent")
cowplot::plot_grid(v1)
dev.off()











