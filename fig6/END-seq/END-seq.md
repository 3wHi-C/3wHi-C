## 1. END-seq Data Processing Script Usage Guide

Script Overview This is a Bash script designed for processing END-seq data, which includes the following steps:

1.  Quality control and adapter trimming (using Trim Galore)
2.  Sequence alignment (using Bowtie2)
3.  BAM file sorting (using Samtools)
4.  Duplicate marking and removal (using Picard)

## 2. Usage Instructions

bash ENDseq_processing.sh <working_directory_path> <sample1> <sample2> [sample3...]
