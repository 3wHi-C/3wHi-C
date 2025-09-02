library(data.table)
library(foreach)
library(doParallel)

setwd("/mnt/disk1/1/gy/Cwalks/new_shuf")
dt <- fread("k562.cadj.annot.txt", col.names = c("chr1", "pos1", "chr2", "pos2", "category"), header = TRUE)

grouped_data <- split(dt, dt$category)
registerDoParallel(cores = 40)

# 并行处理每个分组
data_final <- foreach(
  i = 1:length(grouped_data),
  .combine = rbind,
  .packages = c("data.table")
) %dopar% {
  group <- grouped_data[[i]]
  if (nrow(group) < 3) return(NULL)
  group[, from := paste(chr1, pos1, sep = ":")]
  group[, to := paste(chr2, pos2, sep = ":")]
  num_stats <- table(c(group$from, group$to))
  num_stats <- num_stats[num_stats == 1]
  if (length(num_stats) < 2) return(NULL)
  start_pos <- names(num_stats)[1]
  end_pos <- names(num_stats)[2]
  connection_map <- setNames(group$to, group$from)
  current_point <- start_pos
  path <- current_point
  while (current_point %in% names(connection_map)) {
    next_point <- connection_map[current_point]
    path <- c(path, next_point)
    current_point <- next_point
  }

  if (length(path) >= 2) {
    data.frame(V1 = path, group = names(grouped_data)[i], stringsAsFactors = FALSE)
  } else {
    NULL
  }
}

stopImplicitCluster()

write.table(data_final, file = "/mnt/disk1/1/gy/Cwalks/new_shuf/zheng.pos.txt2", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

