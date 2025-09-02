#!/bin/bash

cd /mnt/disk2/3/gy/3w_fig3/MD/heteropolymers_meng_PNAS/50kb_MD_binomial_result/single_2w_3w

for i in {1..500}
do
	awk -v OFS="\t" '{if($1=="cell_"a)print $2-1","$3-1"\n"$3-1","$2-1}' a=$i ../chr1_pairs_distances_500models_lt1.5.txt2 > cliques.${i}.cells
	line_count=$(wc -l < cliques.${i}.cells)
	max_value=$(awk -F ',' '{print $1"\n"$2}' cliques.${i}.cells | sort -n | tail -n 1)
	echo "$((max_value + 1))" > tmp_file
	echo "$line_count" >> tmp_file
	cat cliques.${i}.cells >> tmp_file
	/mnt/disk1/1/gy/software/quick-cliques/bin/qc --input-file=tmp_file --algorithm=adjlist > final_cell_${i}_cliques
	cat final_cell_${i}_cliques | awk -F " " '{if(NF==2)print $0}' | sed '1,2d' > final_only_2w_cell_${i}_cliques
	cat final_cell_${i}_cliques | awk -F " " '{if(NF==3)print $0}' | sed '1,2d' > final_only_3w_cell_${i}_cliques
	cat final_only_2w_cell_${i}_cliques | awk -F " " -v a=${i} '{print "cell_"a" "$0}' >> total_only_2w_500cell_lt1.5
	cat final_only_3w_cell_${i}_cliques | awk -F " " -v a=${i} '{print "cell_"a" "$0}' >> total_only_3w_500cell_lt1.5
	rm -rf final_cell_* cliques.${i}.cells final_only_2w_cell_${i}_cliques final_only_3w_cell_${i}_cliques
done
