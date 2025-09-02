##RegionR
library("regioneR")

filelist <- c("/mnt/disk6/3/3w_manuscript/sub_boundary/method2_call_boundary/Total/3w_s_new/permutation_test/cell/CTCF_boundary/3w_s_CTCF_boundary.txt",
              "/mnt/disk6/3/3w_manuscript/sub_boundary/method2_call_boundary/Total/3w_s_new/permutation_test/cell/CTCF_boundary/Common_CTCF_boundary.txt",
              "/mnt/disk6/3/3w_manuscript/sub_boundary/method2_call_boundary/Total/3w_s_new/permutation_test/cell/CTCF_boundary/2w_s_CTCF_boundary.txt")

cell_list <- paste("/mnt/disk2/3/gy/cell_type_hic/X_chip_data/TAD_10kb_arrowhead_boundaries.",1:10,"cell.bed",sep = "")

result_tab <- data.frame()
for (i in 1:length(filelist)) {
  data_file <- read.table(filelist[i],header = F)
  query_regions <- toGRanges(data_file)
  for (t in 1:length(cell_list)) {
    cell_type <- read.table(cell_list[t],header = F)
    target_regions <- toGRanges(cell_type[,1:3])
    set.seed(123)
    result <- permTest(
      A = query_regions,
      B = target_regions,
      ntimes = 500,
      genome = "hg38",
      evaluate.function = numOverlaps,
      randomize.function = randomizeRegions,
      force.parallel = 40,
      mask = getGenomeAndMask("hg38"))
    result_tmp <- c(i,t,result$numOverlaps$zscore,result$numOverlaps$pval,result$numOverlaps$observed)
    result_tab <- rbind(result_tab,result_tmp)
  }
}
?permTest
colnames(result_tab) <- c("Type","cell","zscore","P","O")
write.table(result_tab,
file = "/mnt/disk6/3/3w_manuscript/sub_boundary/method2_call_boundary/Total/3w_s_new/permutation_test/cell/CTCF_boundary/permutation_test_result.txt",
col.names = T,row.names = F,quote = F,sep = "\t")




