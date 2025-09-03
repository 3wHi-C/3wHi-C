#!/bin/bash

computeMatrix=/opt/deepTools-3.3.1/bin/computeMatrix
#plot=/opt/deepTools-3.3.1/bin/plotProfile
plot=/opt/deepTools-3.3.1/bin/plotHeatmap
WORKDIR=/mnt/disk2/3/gy/fig6/change
cd $WORKDIR

$computeMatrix reference-point -S /mnt/disk1/1/gy/hhy/END-seq/merge/break_point_bw/total_APH_DSB.bw /mnt/disk1/1/gy/hhy/END-seq/merge/break_point_bw/total_DMSO_DSB.bw \
	-R /mnt/disk1/1/gy/hhy/END-seq/merge/peak/broad/APH_specific_broad_summit.bed \
	--referencePoint center \
	--beforeRegionStartLength 5000 \
         --afterRegionStartLength 5000 \
	--binSize 100 \
	--missingDataAsZero \
	--numberOfProcessors 50 \
	--samplesLabel APH DMSO \
	--outFileName APH_s_end_seq.gz \
	--outFileSortedRegions 1.bed

$plot --matrixFile APH_s_end_seq.gz \
	--outFileName APH_s_end_seq.pdf \
	--refPointLabel "Peak summit" \
	--dpi 400 \
	--samplesLabel APH DMSO \
	--dpi 400 \
	--colorMap Blues \
	--heatmapHeight 8 \
	--zMin 0 --zMax 0.003 \
	--yMin 0 --yMax 0.002

