##cells cliques
setwd("/mnt/disk2/3/gy/3w_fig3/MD/heteropolymers_meng_PNAS")

library("doParallel")
library("foreach")
library("dplyr")
library("data.table")

cell_data <- fread("chr1_pairs_distances_500models_lt2.5.txt2",header = F)
cell_data$V2 <- as.numeric(cell_data$V2)
cell_data$V3 <- as.numeric(cell_data$V3)
##cliques
cell_list <- unique(cell_data$V1)

n_cores <- 15
cl <- makeCluster(n_cores)
registerDoParallel(cl)

foreach(current_cell = cell_list, 
        .packages = c("dplyr"),
        .export = "cell_data") %dopar% {
          tmp_cell_data <- cell_data[cell_data$V1 == current_cell,]
          result <- tmp_cell_data %>% 
            inner_join(tmp_cell_data[,2:3],by=c("V2"="V2"),relationship = "many-to-many") %>%
            inner_join(tmp_cell_data[,2:3],by=c("V3.x"="V2","V3.y"="V3"))
          
          output_file <- paste0("/mnt/disk2/3/gy/3w_fig3/MD/heteropolymers_meng_PNAS/matrix3/",current_cell, "_matrix3.txt")
          write.table(result[,c(1,2,3,5)],file = output_file,sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
        }

stopCluster(cl)


