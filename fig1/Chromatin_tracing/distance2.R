data <- read.table("K562_chr21-28-30Mb.csv2",header = F,sep = ",")
data <- na.omit(data)
library(parallel)
max_cores <- 30

compute_distance_parallel <- function(group_data) {
  result <- data.frame()
  
  for (i in 1:(nrow(group_data) - 1)) {
    for (j in (i + 1):nrow(group_data)) {
      probe1 <- group_data[i,]
      probe2 <- group_data[j,]
      distance <- sqrt((probe1$V3 - probe2$V3)^2 + 
                         (probe1$V4 - probe2$V4)^2 + 
                         (probe1$V5 - probe2$V5)^2)     
      result <- rbind(result, data.frame(
        group = group_data$V1[1],
        probe_pair = paste(probe1$V2, probe2$V2, sep = "-"),
        euclidean_distance = distance
      ))
    }
  }
  return(result)
}

compute_group_distances <- function(df, max_cores) {
  groups <- missdata
  
  group_results <- mclapply(groups, function(group) {
    group_data <- df[df$V1 == group, ]
    compute_distance_parallel(group_data)
  }, mc.cores = max_cores)
  
  result <- do.call(rbind, group_results)
  
  return(result)
}

distance_result <- compute_group_distances(raw_data,max_cores)
write.table(distance_result, "/mnt/disk6/3/software/ChromatinImaging-master/Data/probe_distances_parallel.csv2",row.names = F, col.names = F, quote = F, sep = "\t")


