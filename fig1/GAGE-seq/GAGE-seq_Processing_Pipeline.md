# GAGE-seq Data Processing Pipeline

This pipeline processes raw GAGE-seq data, filters interactions, and performs clique analysis. Adapted from K562 single-cell Hi-C data (GEO accession: GSM7657692).

## Dependencies

-   Unix/Linux environment
-   Bash shell
-   GNU Core Utilities
-   R (for binning steps)
-   `cliques.sh` (custom clique analysis tool)

------------------------------------------------------------------------

## Pipeline Steps

### 1. Raw Data Preparation

\`\`\`bash \# Decompress the raw pairs file gunzip GSM7657692_contact_K3-0208-PL2-HiC.pairs.gz

# Remove non-chromosomal entries and format columns

cat GSM7657692_contact_K3-0208-PL2-HiC.pairs \| grep -v "mm10" \| awk -vOFS="\t" '{print substr(\$1,6),\$2,substr(\$3,6),\$4,\$5}' \> K562_GSM7657692_contact_K3-0208-PL2-HiC.pairs

# Generate cell interaction statistics

cat K562_GSM7657692_contact_K3-0208-PL2-HiC.pairs \| cut -f5 \| sort -k1,1 -S 30% --parallel=40 \| uniq -c \> cell_stats

# Retain only cells with ≥10,000 interactions (using well.ID as filter)

grep -w -f well.ID K562_GSM7657692_contact_K3-0208-PL2-HiC.pairs \> filter_K562_GSM7657692_contact_K3-0208-PL2-HiC.pairs

# Filter for cis-chromosomal interactions

cat filter_K562_GSM7657692_contact_K3-0208-PL2-HiC.pairs \| awk -v OFS="\t" '{if(\$1==\$3)print \$0}' \> cis_filter_K562_GSM7657692_contact_K3-0208-PL2-HiC.pairs cat

cis_filter_K562_GSM7657692_contact_K3-0208-PL2-HiC.pairs \| awk -v OFS="\t" '{if(\$4\>\$2)print \$0;else print \$1,\$4,\$3,\$2,\$5}' \> cis_filter_K562_GSM7657692_contact_K3-0208-PL2-HiC.pairs2

# Bin interactions using R script (GAGE_process.R)

# Converts pairwise interactions to binned matrix format

# Perform maximal clique analysis (requires cliques.sh)

# Output: total_cell_cliques

# Filter for three-way cliques

cat total_cell_cliques \| awk -v OFS="\t" '{if(NF==4)print \$0}' \> total_cell_cliques_3w
