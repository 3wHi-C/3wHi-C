#!/bin/bash
if [ $# -lt 2 ]; then
    echo "Usage: $0 <directory> <sample1> <sample2> ..."
    exit 1
fi

trim_galore=/usr/local/bin/trim_galore
cutadapt=/usr/local/bin/cutadapt
bowtie2=/opt/bowtie2-2.4.5/bowtie2
JAVA=/usr/bin/java
PICARD=/mnt/disk1/1/gy/software/Picard/picard.jar
samtools=/usr/local/bin/samtools

WORK_DIR="$1"
shift
SAMPLES=("$@")

cd "$WORK_DIR" || exit 1

mkdir -p trimAdapt

for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: $sample"
    
    $trim_galore -q 20 -j 10 --phred33 --stringency 3 \
        --length 20 -e 0.1 --fastqc --gzip \
        -o trimAdapt/ \
        --path_to_cutadapt $cutadapt \
        "${sample}_R1.fq.gz"
    
    $bowtie2 -t -q -p 60 --very-sensitive -X 2000 \
        --score-min L,0,-0.4 \
        -x /ssd/index/bowtie2/hg38XX \
        -U "trimAdapt/${sample}_R1_trimmed.fq.gz" \
        2> "${sample}.stat" | \
        $samtools view -@ 60 -bh -q 10 -o "${sample}_q10.bam" -
    
    $samtools sort -@ 40 -o "${sample}_q10.sort.bam" "${sample}_q10.bam"
    rm "${sample}_q10.bam"
    
    $JAVA -Xmx40g -jar $PICARD MarkDuplicates \
        I="${sample}_q10.sort.bam" \
        O="${sample}_dedu.bam" \
        M="${sample}.m" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true
done

echo "All samples processed successfully!"
