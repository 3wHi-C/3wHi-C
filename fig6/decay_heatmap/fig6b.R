setwd("/mnt/disk2/3/gy/fig6/expected")
APH <- fread("/mnt/disk2/3/gy/APH_DMSO/APH/100kb_m3/Total_APH_intra_rawcount.matrix3",header = F)
DMSO <- fread("/mnt/disk2/3/gy/APH_DMSO/DMSO/100kb_m3/Total_DMSO_intra_rawcount.matrix3",header = F)

DMSO$distance1 <- DMSO$V2 - DMSO$V1
DMSO$distance2 <- DMSO$V3 - DMSO$V2

expected_DMSO <- DMSO %>% group_by(distance1,distance2) %>% summarise(sum=sum(V4))
expected_DMSO$P <- expected_DMSO$sum/sum(expected_DMSO$sum)

expected_DMSO$distance1 <- expected_DMSO$distance1*100000 + 1
expected_DMSO$distance2 <- expected_DMSO$distance2*100000 + 1

expected_DMSO$bin_x <- cut(log10(expected_DMSO$distance1),breaks = c(0,seq(5.0000,8.5,by=0.5)), include.lowest = TRUE)
expected_DMSO$bin_y <- cut(log10(expected_DMSO$distance2),breaks = c(0,seq(5.0000,8.5,by=0.5)), include.lowest = TRUE)

##APH
stats_APH <- expected_APH[,c(5,6,4)] %>% group_by(bin_x,bin_y) %>% summarise(sum=sum(P))

##DMSO
stats_DMSO <- expected_DMSO[,c(5,6,4)] %>% group_by(bin_x,bin_y) %>% summarise(sum=sum(P))

stats_APH_DMSO <- left_join(stats_DMSO,stats_APH,by=c("bin_x"="bin_x","bin_y"="bin_y"))
colnames(stats_APH_DMSO)[3:4] <- c("DMSO","APH")
stats_APH_DMSO$delta <- stats_APH_DMSO$APH - stats_APH_DMSO$DMSO
#stats_APH_DMSO$delta <- log2(stats_APH_DMSO$APH/stats_APH_DMSO$DMSO)

pdf("/mnt/disk2/3/gy/fig6/expected/APH_DMSO_delta_decay_100kb.pdf", width = 4.5, height = 3.2, bg = "transparent")
ggplot(data = stats_APH_DMSO,mapping = aes(x = bin_x, y = bin_y))+
  theme_classic()+
  #geom_tile(aes(fill=log10(abs(delta))), colour=NA)+
  #geom_tile(aes(fill=log10(Expected)), colour=NA)+
  geom_tile(aes(fill=delta), colour=NA)+
  coord_equal()+
  scale_fill_gradient2(low = "blue",mid = "white", high = "red",name="Delta expected",midpoint = 0,limits = c(-0.1,0.1))+
  #scale_fill_gradient(low = "white", high = "red",name="log10\n(expected)",limits = c(-8,0))+ #limits = c(-8,0)
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15,angle = 90,colour = "black"),
        axis.text.y = element_text(size = 15,colour = "black"),
        #legend.position = c(0.7, 0.3),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 15, colour = 'black'),
        legend.text = element_text(size = 15, colour = 'black'))+
  xlab("")+ylab("")
dev.off()
