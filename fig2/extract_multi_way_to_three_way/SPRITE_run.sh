#!/bin/bash

cd /home/gy/SPRITE/GSE114242_human_combined_clusters

awk -F '[:\t]' '{delete count; for(i=2; i<=NF; i+=2) count[$i]++; printf "%s",$1; for(type in count) printf "\t%s:%d", type, count[type]; printf "\n"}' human.combined.mapq-ge30.clusters > line_count.txt
awk -F '[:\t]' '{for(i=3; i<=NF; i+=2){ if($i >= 3 && $i <=10 && i >= 3) printf "%s\t%s\t%s\n",$1,$(i-1),$i}}' line_count.txt > line_count_3_10.txt
awk 'BEGIN{FS=OFS="\t"} NR==FNR{reads[$1]=$2;next} $1 in reads{split($0,array," "); for(i=2;i<=length(array);i++){split(array[i],chr_pos,":");if(chr_pos[1]==reads[$1]){print $1, array[i]}}}' chr21_line_count_3_10.txt human.combined.mapq-ge30.clusters > chr21_select_3_10_human.combined.mapq-ge30.clusters
