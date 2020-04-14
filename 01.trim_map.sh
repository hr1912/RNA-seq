#! /bin/bash

java='/usr/bin/java'
trim='/path/trimmomatic-0.38.jar'
adapter='/path/illumina/Adapter_PE.fa'

sample=$1
fq1=$2
fq2=$3
filt_fn_r1="${sample}_trimmed.R1.fastq.gz"
filt_fn_r2="${sample}_trimmed.R2.fastq.gz"
unp_fn_r1="${sample}_dumped.R1.fastq.gz"
unp_fn_r2="${sample}_dumped.R2.fastq.gz"

dir_out='.'
mkdir ${dir_out}/${sample}

${java} -jar ${trim} PE ${fq1} ${fq2} \
		${dir_out}/${sample}/${filt_fn_r1} \
		${dir_out}/${sample}/${unp_fn_r1} \
		${dir_out}/${sample}/${filt_fn_r2} \
		${dir_out}/${sample}/${unp_fn_r2} \
        'ILLUMINACLIP:'${adapter}':2:30:10' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 2> ${dir_out}/${sample}/read_qc_stat.txt

idx_human='/path/index/gencode_v28'
anno_human='/path/index/gencode.v28.annotation.gtf'

hisat2 -p 4 --dta -x $idx_human \
	-1 ${dir_out}/${sample}/$filt_fn_r1 \
	-2 ${dir_out}/${sample}/$filt_fn_r2 \
   | samtools view -bS - > ${dir_out}/${sample}/${sample}.bam

samtools sort ${dir_out}/${sample}/${sample}.bam "${sample}.sorted"

mv "${sample}.sorted.bam" ${dir_out}/${sample}

