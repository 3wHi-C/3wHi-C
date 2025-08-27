#!/bin/bash
get_script_dir() {
    echo "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
}
script_dir="$(get_script_dir_simple)"
MAX_CPU_NUMBER=`cat /proc/cpuinfo |grep "processor" | wc -l`
processingOneInstance(){
	prefix=${1}
	if [ ${3} ]; then
		v=1
	else
		v=0
	fi
	cd ${1}
	tmp_fifofile='${1}.fifo'
	mkfifo  $tmp_fifofile
	exec 7<>$tmp_fifofile
	rm  $tmp_fifofile
	for (( i=0 ;i < `expr $MAX_CPU_NUMBER / ${2}` ;i++ )); do
		echo
	done >&7
	if [ ! -f config* ];then
		cp ${script_dir}/config-general.txt config-${1}.txt
	fi
	configFile=`ls config* | head -n 1`
	if [ ! -f UMI.map ];then
		cp ${script_dir}/UMI.map ../UMI.map
	fi
	mkdir blast split_fa
	fastp -i ${1}_R1.fq.gz -I ${1}_R2.fq.gz -m --merged_out merged_${1}.fq.gz
	${script_dir}/fq_to_fa merged_${1}.fq.gz
	split -l 200000 -d merged_${1}.fa split_fa/x
	for i in `ls split_fa`; do
		read -u7
		{
			blastn -query split_fa/${i} -subject ${script_dir}/3w3e.fa -word_size 20 -outfmt 6 -out blast/${i}.tmp1 
			echo  >&7
		}&
	done
	wait
	exec  7>&-
	
	#for i in {0..9};do {
	#	for j in `ls split_fa/*${i}`;do {
	#		blastn -query ${j} -subject ${script_dir}/3w3e.fa -word_size 30 -outfmt 6 -out $blast/${j##*/}.tmp1
	#	}; done
	#}& done
	#wait
	cat blast/*.tmp1 > blast/${1}.blast
	rm -f blast/*.tmp1
	rm -rf split_fa
	python3 ${script_dir}/mark_strand.py ${1}
	python3 ${script_dir}/DNA_binarization.py ${1}
	python3 ${script_dir}/split_index.py ${1}
	cd hicpro
	hic-pro -s mapping -s proc_hic -s quality_checks -c ../config-${1}.txt -i rawdata -o ${1}_hicpro
	cd ..
	samtools view hicpro/${1}_hicpro/hic_results/data/${1}/${1}_hg38XX.bwt2pairs_interaction.bam > ${1}_hg38XX.bwt2pairs_interaction.sam
	shift 3
	while [ $# -gt 0 ]; do
		samtools view ${1}/${prefix}/hicpro/${prefix}_hicpro/hic_results/data/${prefix}/${prefix}_hg38XX.bwt2pairs_interaction.bam >> ${prefix}_hg38XX.bwt2pairs_interaction.sam
		cat ${1}/${prefix}/${prefix}_R1.umi >> ${prefix}_R1.umi
		cat ${1}/${prefix}/${prefix}_R2.umi >> ${prefix}_R2.umi
		cat ${1}/${prefix}/${prefix}.1.l.marker >> ${prefix}.1.l.marker
		shift 1
	done
	python3 ${script_dir}/sam_to_w3.py ${prefix} ${v}
}

show_help() {
    cat << EOF
Usage: $SCRIPT_NAME [options]


Options:
  -h, --help          Show help message;
  -s, --single        Run script in a directory containing fastq files of single library;
  -p, --threads INT   Number of parallel samples;
  -f, --force         Clean directory before analysis;
  -r, --rep <PATH>    Processed technical duplicates directory. 

示例:
  $SCRIPT_NAME -f input.txt -o ./output
  $SCRIPT_NAME --help


EOF
}

test=0
p=1
s=''
rep=()

if [ $# -eq 0 ]; then
    show_help
    exit 0
fi

while [ $# -gt 0 ]; do
	case ${1} in
	-h|--help)
    show_help
    exit 0
	-s|--single)
		shift 1
		s=${1}
		shift 1
		;;
	-p|--threads)
		shift 1
		p=${1}
		shift 1
		;;
	-f|--force)
		rm -rf */blast/ */config-*  */hicpro/ */merged* */fastp* */split_fa
		shift 1
		;;
	-r|--rep)
		shift 1
		if [ ${1:0:1} != '-' ];then
			rep+=`realpath ${1}`
			shift 1
		fi
		;;
	*)
		shift
		;;
	esac
done

if [ '$s' = '' ]; then
	processingOneInstance $s $p $test
else
	tmp_fifofile='/tmp/$$.fifo'
	mkfifo $tmp_fifofile
	exec 6<>$tmp_fifofile
	rm  $tmp_fifofile
	for (( i=0 ;i < $p ;i++ )); do
		echo
	done >&6
	for i in `ls`; do
		if [ -d  ${i} ]; then
			read  -u6
			{
				processingOneInstance ${i} $p $test ${rep[*]}
				echo  >&6
			}&
		fi
	done
	wait
	exec  6>&-
fi
wait

${script_dir}/financy_curve.py
