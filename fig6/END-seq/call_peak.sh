#!/bin/bash

OUTPUTDIR=/mnt/disk1/1/gy/hhy/END-seq/merge/peak/broad
INPUTDIR=/mnt/disk1/1/gy/hhy/END-seq/merge

MACS2=/opt/MACS-2.2.5/bin/macs2

for i in total_APH_sort total_DMSO_sort
do
	$MACS2 callpeak -t $INPUTDIR/${i}.bam \
		--name ${i} \
		-B \
		-f BAM \
		-g hs \
		--keep-dup all \
		--outdir $OUTPUTDIR \
		-q 0.05 \
		--nolambda \
		--nomodel \
		--broad
done
