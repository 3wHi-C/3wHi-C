setwd("/mnt/disk6/3/software/ChromatinImaging-master/Data/mydata")
library("dplyr")
data <- readLines("/mnt/disk6/3/software/ChromatinImaging-master/Data/total_cell_cliques")
split_strings <- strsplit(data," ")
#data_final <- sapply(split_strings, function(x) paste(x[-1], collapse = " "))
data_final <- sapply(split_strings, function(x) {
  sorted_numbers <- sort(as.numeric(x[-1]))
  paste(sorted_numbers, collapse = " ")
})

#count_matches_1 <- sapply(prefixes, function(prefix) sum(grepl(prefix,data_final)))
count_matches <- sapply(data_final, function(prefix) {
  numbers <- unlist(strsplit(prefix, " "))
  regex <- paste(sapply(numbers, function(num) paste0("\\b",num,"\\b")), collapse = ".*")
  sum(grepl(regex,data_final))
})

aa <- as.data.frame(data_final)
aa$prefix_count <- count_matches

write.table(aa,file = "complex_in_image.results",sep = ",",quote = F,row.names = F,col.names = F)
#write.table(aa,file = "complex_in_image.results",sep = "\t",quote = F,row.names = F,col.names = F)