#! /bin/bash
PLASMID_FINDER='/home/cxz163430/software/plasmidfinder'
ASM='/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT'
HIP='/ms/11/cong/project/HiC/data/ref/HIP11704'
V583='/ms/11/cong/project/HiC/data/ref/V583'
PLS='/ms/11/cong/data/PLS_db/1903/2019_03_05.fna'

file='/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/coassembly'
BASE='Coassembly'

blastn -query ${file}/final.contigs.fa -db ${V583}/V583 -perc_identity 90 -num_threads 8 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/${BASE}_V583_genome.tsv
blastn -query ${file}/final.contigs.fa -db ${HIP}/HIP11704 -perc_identity 90 -num_threads 8 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/${BASE}_HIP_genome.tsv
blastn -query ${file}/final.contigs.fa -db ${PLASMID_FINDER}/plasmidfinder_db/AllRep -perc_identity 90 -num_threads 8 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/${BASE}_rep.tsv
blastn -query ${file}/final.contigs.fa -db ${PLS} -perc_identity 90 -num_threads 8 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/${BASE}_plsdb.tsv

awk '$8 >=800' ${BASE}_V583_genome.tsv > ${BASE}_V583_90_800.tsv
awk '$8 >=800' ${BASE}_HIP_genome.tsv > ${BASE}_HIP_90_800.tsv
awk '$8 >=800' ${BASE}_rep.tsv > ${BASE}_rep_90_800.tsv
