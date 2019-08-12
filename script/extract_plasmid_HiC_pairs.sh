#! /bin/sh
##use .script plasmid_sequence HiC_1 HiC_2 output name

HIP_plasmid=$1
HiC_1=$2
HiC_2=$3
output=$4

bwa index  ${HIP_plasmid}
bwa mem -t 24 -5SP ${HIP_plasmid} ${HiC_1} ${HiC_2} | 
#    samblaster | 
	samtools view -h -F 2316 -bS -@ 24 - | \
        samtools sort -n -o ${output}_2316.bam -
samtools bam2fq ${output}_2316.bam | seqtk seq -A > ${output}_2316.fa


bwa mem -t 24 -5SP ${HIP_plasmid} ${HiC_1} ${HiC_2} |
#    samblaster | 
	samtools view -h -f 8 -F 2308 -bS -@ 24 - | \
	samtools sort -n -o ${output}_2308.bam -
	samtools bam2fq ${output}_2308.bam | seqtk seq -A > ${output}_2308.fa

bwa mem -t 24 -5SP ${HIP_plasmid} ${HiC_1} ${HiC_2} |
#    samblaster | 
        samtools view -h -f 4 -F 2312 -bS -@ 24 - | \
        samtools sort -n -o ${output}_2312.bam -
samtools bam2fq ${output}_2312.bam | seqtk seq -A > ${output}_2312.fa
