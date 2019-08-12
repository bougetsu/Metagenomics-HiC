#! /bin/bash

PLASMID_FINDER='/home/cxz163430/software/plasmidfinder'

#makeblastdb -in ${PLASMID_FINDER}/plasmidfinder_db/gram_positive.fsa -dbtype nucl -parse_seqids -out ${PLASMID_FINDER}/plasmidfinder_db/GramPos
HIP='/ms/11/cong/project/HiC/data/ref/HIP11704'
OG1='/ms/11/cong/project/HiC/data/ref/OG1RF/GCF_004006275.1_ASM400627v1_genomic.fna'
T11='/ms/11/cong/project/HiC/data/ref/T11/GCF_000157475.1_ASM15747v1_genomic.fna'
PLS='/ms/11/cong/data/PLS_db/1903/2019_03_05.fna'
contig='/ms/11/cong/project/HiC/data/raw_data/1st_mouse/ASSEMBLY/IDBAUD/contig.fa'
RES='/ms/11/cong/project/HiC/data/raw_data/1st_mouse/blast_result'
pAD1='/ms/11/cong/project/HiC/data/ref/plasmid/pAD1.fasta'

blastn -query ${contig} -db ${OG1} -perc_identity 90 -num_threads 8  -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${RES}/idbuad_OG1.tsv
blastn -query ${contig} -db ${T11} -perc_identity 90 -num_threads 8  -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${RES}/idbuad_T11.tsv
blastn -query ${contig} -db ${pAD1} -perc_identity 90 -num_threads 8  -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${RES}/idbuad_pAD1.tsv
blastn -query ${contig} -db ${PLASMID_FINDER}/plasmidfinder_db/AllRep -perc_identity 90 -num_threads 8  -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${RES}/idbuad_rep.tsv
blastn -query ${contig} -db ${PLS} -perc_identity 90 -num_threads 8  -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${RES}/idbuad_pls.tsv