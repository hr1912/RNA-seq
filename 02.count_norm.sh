sample=$1

py27='/usr/bin/python'
dir_out='.'

if [-d ${dir_out}/${sample}] 
then
	echo "Directory ${dir_out}/${sample} exists."
else
	mkdir ${dir_out}/${sample}
fi

idx_human='/path/index/gencode_v28'
anno_human='/path/index/gencode.v28.annotation.gtf'


stringtie ${sample}.sorted.bam \
	-e -B -p 4 -G ${anno_human} \
	-o ${dir_out}/${sample}/${sample}.table.txt \
	-A ${dir_out}/${sample}/${sample}.genelevel_table.txt
