#!/bin/bash

RepeatMasker=/mnt/disk1/1/gy/software/RepeatMasker/RepeatMasker/RepeatMasker

cd /mnt/disk6/3/useful_file/Transposable_elements/RepeatMasker_my

${RepeatMasker} -pa 20 -s -nolow -species human -gccalc -gff -cutoff 200 -no_is -dir ./RepeatMasker_out /mnt/disk1/1/gy/useful-file/genome/hg38F.fna

