# Multi-way Interaction Extraction and Annotation Pipeline

This repository contains scripts to extract pairwise/three-way interactions from multi-way data and annotate them with radial layer information.

## 1. Extract Three-way Interactions

\`\`\`bash python3 extract_three_way.py

## 2. Data Sources

## HiPore-C Data

• ​​Location​​: /mnt/disk1/1/gy/HiPore-C/distance • ​​Example​​: Chromosome 1 data • ​​Format​​: Directly usable after extraction

## C-walks Data

• ​​Location​​: /mnt/disk1/1/gy/Cwalks/final_new_shuf • ​​Preprocessing Required​​: Run adjust_Cwalks_reads.R to rearrange fragment order Convert coordinates from hg19 to hg38

## SPRITE Data

• ​​Location​​: /mnt/disk1/1/gy/SPRITE/GSE114242_human_combined_clusters/hg38_final_3w • ​​Processing Pipeline​​: bash SPRITE_run.sh

## This script:

a.  Counts chromosome complexes per read (3-10 range)
b.  Formats output as reads chr:pos
c.  Randomly selects three regions per read

## 3. Radial Layer Annotation

bedtools intersect -a \<fragments.bed\> -b \<radial_layers.bed\> -wa -wb \> annotated.bed
