##Perform Spontaneous model analysis and binomial test for each chromosome
##====================================load package=============================================
library("data.table")
library("dplyr")
library("stats")
library("parallel")
##====================================prepare data==============================================
##res: 50Kb 75Kb 100Kb
bed_50kb <- fread("/mnt/disk6/3/useful_file/K562_microC/matrix/K562_MicroC_R1+R2+R3+R4.75000_abs.bed",header = F)
bed_50kb_stats <- bed_50kb %>% group_by(V1) %>% summarise(min=min(V4),max=max(V4))
chr_list <- paste("chr",c(1:22,"X"),sep = "")
#chr_list <- c("chr1","chr2")
matrix_50kb <- fread("/mnt/disk6/3/useful_file/K562_microC/matrix/K562_MicroC_R1+R2+R3+R4.75000_iced.matrix",header = F)
matrix3_50kb <- fread("/mnt/disk6/3/3w_manuscript/all_m3/75kb_matrix3/nor_B_cis_gt500_75kb.matrix3",header = F)
matrix3_50kb$V4 <- round(matrix3_50kb$V4)
for (chr in 1:length(chr_list)) {
  start <- bed_50kb_stats[bed_50kb_stats$V1==chr_list[chr],]$min
  end <- bed_50kb_stats[bed_50kb_stats$V1==chr_list[chr],]$max
  matrix_50kb_tmp <- matrix_50kb[matrix_50kb$V1>=start & matrix_50kb$V1<=end & matrix_50kb$V2>=start & matrix_50kb$V2<=end,]
  matrix_50kb_tmp$w2_probe <- matrix_50kb_tmp$V3/sum(matrix_50kb_tmp$V3)
  matrix3_50kb_tmp <- matrix3_50kb[matrix3_50kb$V1>=start & matrix3_50kb$V1<=end & 
                                     matrix3_50kb$V2>=start & matrix3_50kb$V2<=end & 
                                     matrix3_50kb$V3>=start & matrix3_50kb$V3<=end,]
  matrix3_50kb_tmp <- inner_join(matrix3_50kb_tmp,matrix_50kb_tmp[,c(1,2,4)],by=c("V1"="V1","V2"="V2"))
  matrix3_50kb_tmp <- inner_join(matrix3_50kb_tmp,matrix_50kb_tmp[,c(1,2,4)],by=c("V1"="V1","V3"="V2"))
  matrix3_50kb_tmp <- inner_join(matrix3_50kb_tmp,matrix_50kb_tmp[,c(1,2,4)],by=c("V2"="V1","V3"="V2"))
  expected_probe <- sapply(1:nrow(matrix3_50kb_tmp), function(i) {
    x <- matrix3_50kb_tmp$w2_probe.x[i]
    y <- matrix3_50kb_tmp$w2_probe.y[i]
    z <- matrix3_50kb_tmp$w2_probe[i]
    geometric_mean <- sqrt(x * y * z)
    part1 <- geometric_mean
    part2 <- min(x, y, z)
    min(part1, part2)
  })
  matrix3_50kb_tmp$expected_probe <- expected_probe
  matrix3_50kb_tmp$expected_probe <- matrix3_50kb_tmp$expected_probe/sum(matrix3_50kb_tmp$expected_probe)
  matrix3_50kb_tmp$sampleN <- sum(matrix3_50kb_tmp$V4)
  matrix3_50kb_tmp$expected_count <- matrix3_50kb_tmp$expected_probe*matrix3_50kb_tmp$sampleN
  output_name <- paste("/mnt/disk2/3/gy/3w_fig3/geometric_mean/75kb_nor_binomial_result/",chr_list[chr],"_matrix3_75kb_expected.txt",sep = "")
  write.table(matrix3_50kb_tmp,file = output_name,sep = "\t",col.names = F,row.names = F,quote = F)
  

  total_rows <- nrow(matrix3_50kb_tmp)
  chunk_size <- 50000
  num_chunks <- ceiling(total_rows/chunk_size)
  for (t in 1:num_chunks) {
    start_row <- (t-1)*chunk_size+1
    end_row <- min(t*chunk_size,total_rows)
    matrix_10kb_test <- matrix3_50kb_tmp[start_row:end_row,]
    p_values <- mclapply(seq(nrow(matrix_10kb_test)), function(i) {
      row <- matrix_10kb_test[i,]
      tryCatch({
        binom.test(as.numeric(row$V4), as.numeric(row$sampleN),as.numeric(row$expected_probe),
                   alternative = "two.sided")$p.value
      }, error = function(e) NA)
    }, mc.cores = 10)
    p_values <- unlist(p_values)
    matrix_10kb_test$P_value <- p_values
    output_bio_file <- paste("/mnt/disk2/3/gy/3w_fig3/geometric_mean/75kb_nor_binomial_result/",chr_list[chr],"_w3_P_expected_",t,".txt",sep = "")
    write.table(matrix_10kb_test,file = output_bio_file,col.names = T,row.names = F,quote = F,sep = "\t")
  }
}

