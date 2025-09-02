setwd("/mnt/disk2/3/gy/3w_fig3/MD/heteropolymers_meng_PNAS")
library("data.table")

chr1_10kb_bin <- fread("chr1_10kb_bin.txt",header = F)

repnum <- 500
simulation_data <- fread("Total_chr1_500models.txt",header = F)
simulation_data$bin_ID <- rep(min(chr1_10kb_bin$V4):max(chr1_10kb_bin$V4),repnum)

compute_distance_optimized <- function(df) {
  library("foreach")
  library("doParallel")
  library("data.table")
  cores <- 15
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (!"bin_ID" %in% names(df)) stop("Dataframe must contain 'bin_ID' column")
  
  results <- foreach(group = unique(df$V5), .combine = rbind) %dopar% {
    group_data <- df[df$V5 == group,]
    n <- nrow(group_data)
    if (n < 2) return(NULL)

    dist_matrix <- as.matrix(dist(group_data[,c("V1","V2","V3")]))
    indices <- which(upper.tri(dist_matrix, diag = FALSE) & dist_matrix <= 2.5, arr.ind = TRUE)
    if (nrow(indices) == 0) return(NULL)
    data.frame(
      group = group,
      probe_pair = paste(group_data$bin_ID[indices[,1]], 
                         group_data$bin_ID[indices[,2]], sep = "-"),
      euclidean_distance = dist_matrix[indices]
    )
  }
  
  stopCluster(cl)
  setDF(results)
  return(results)
}

system.time(
  result <- compute_distance_optimized(simulation_data)
)

write.table(result, "/mnt/disk2/3/gy/3w_fig3/MD/heteropolymers_meng_PNAS/chr1_pairs_distances_500models_lt2.5.txt",
            row.names = F,col.names = T,quote = F,sep = "\t")






