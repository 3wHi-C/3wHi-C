## 1. Distance-based Clique Filtering

Calulate Euclidean distance (requires distance2.R)

\`\`\`bash \# Filter cliques where probe distances are \<150nm from probe_distances_parallel.csv2

# 2. Extract unique cell IDs from filtered cliques

cat K562_dis_150nm.cliques \| awk -F ',' '{print \$1}' \| uniq \> cell_ID

# 3. Perform maximal clique analysis (requires cliques.sh)

# 4. Count occurrences of three-way in imaging data

# common_image.R script outputs: complex_in_image.results

# 5. Filter for three-way and remove duplicates

cat complex_in_image.results \| awk -F '[ ,]' -v OFS="\t" '{ if(NF==4) print \$1+1,\$2+1,\$3+1,\$4 }' \| sort -S 30% --parallel=40 \| uniq \> w3_complex_in_image.results
