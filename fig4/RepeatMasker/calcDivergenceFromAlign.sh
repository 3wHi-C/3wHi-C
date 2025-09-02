#!/bin/bash

perl=/usr/bin/perl
calcDivergenceFromAlign=/mnt/disk1/1/gy/software/RepeatMasker/RepeatMasker/util/calcDivergenceFromAlign.pl

cd /mnt/disk6/3/useful_file/Transposable_elements/RepeatMasker_my

${perl} ${calcDivergenceFromAlign} -s /mnt/disk6/3/useful_file/Transposable_elements/RepeatMasker_my/RepeatMasker_out/summary.TE.out /mnt/disk6/3/useful_file/Transposable_elements/RepeatMasker_my/RepeatMasker_out/hg38F.fna.cat.gz

