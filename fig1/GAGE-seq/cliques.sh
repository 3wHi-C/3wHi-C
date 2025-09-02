#!/bin/bash

cd /mnt/disk6/3/useful_file/scHi-C/GAGE-seq

for i in `cat cell_newID`
do
	awk -F ',' '{if($1==a)print $2-1","$3-1}' a=$i K562_GAGE_seq.freq_100kb.txt > cliques.${i}.cells
	line_count=$(wc -l < cliques.${i}.cells)
	max_value=$(awk -F ',' '{print $1"\n"$2}' cliques.${i}.cells | sort -n | tail -n 1)
	echo "$((max_value + 1))" > tmp_file
	echo "$line_count" >> tmp_file
	cat cliques.${i}.cells >> tmp_file
	/mnt/disk1/1/gy/software/quick-cliques/bin/qc --input-file=tmp_file --algorithm=adjlist > final_gt3_${i}_cliques
	cat final_gt3_${i}_cliques | awk -F " " '{if(NF>=3)print $0}' | sed '1,2d' > final_gt3_${i}_cliques2
	cat final_gt3_${i}_cliques2 | awk -F " " -v a=${i} '{print a" "$0}' >> total_cell_cliques_100kb
	rm -rf final_gt3_* cliques.${i}.cells
done
