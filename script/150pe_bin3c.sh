#! /bin/bash

ASM='/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/150PE/'
BAM='/ms/11/cong/project/HiC/BINNING/150PE/bam'
BIN='/ms/11/cong/project/HiC/BINNING/150PE/bin3c'
for file in ${ASM}/*contigs.fa;do
	SAMPLE=${file##*/}
	BASE=${SAMPLE%_150PE.contigs.fa*}
	./bin3C.py mkmap -e MluCI -e Sau3AI -v ${file} ${BAM}/${BASE}_hic2mega2316_150pe.bam ${BIN}/${BASE}
	./bin3C.py cluster --no-spades  -v ${BIN}/${BASE}/contact_map.p.gz ${BIN}/${BASE}_clust
done

