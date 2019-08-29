#! /bin/bash
PLASMID_FINDER='/home/cxz163430/software/plasmidfinder'

#makeblastdb -in ${PLASMID_FINDER}/plasmidfinder_db/gram_positive.fsa -dbtype nucl -parse_seqids -out ${PLASMID_FINDER}/plasmidfinder_db/GramPos
HIP='/ms/11/cong/project/HiC/data/ref/HIP11704'
V583='/ms/11/cong/project/HiC/data/ref/V583'
PLS='/ms/11/cong/data/PLS_db/1903/2019_03_05.fna'
#makeblastdb -in ${V583}/V583.fna -dbtype nucl -parse_seqids -out ${V583}/V583
#makeblastdb -in ${HIP}/HIP11704.fna -dbtype nucl -parse_seqids -out ${HIP}/HIP11704

#blastn -query plasmid_0425.fasta -db ${PLASMID_FINDER}/plasmidfinder_db/GramPosRep -perc_identity 90  -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out test.tsv

#cat ${PLASMID_FINDER}/plasmidfinder_db/gram_positive.fsa ${PLASMID_FINDER}/plasmidfinder_db/enterobacteriaceae.fsa > ${PLASMID_FINDER}/plasmidfinder_db/all_rep.fsa
#makeblastdb -in ${PLASMID_FINDER}/plasmidfinder_db/all_rep.fsa -dbtype nucl -parse_seqids -out ${PLASMID_FINDER}/plasmidfinder_db/AllRep

ASM='/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE'

mkdir -p ${ASM}/blast_result/filtered

file="${ASM}/final.contigs.fa"

#blastn -query ${file} -db ${V583}/V583 -perc_identity 90 -num_threads 16  -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/CoAsm_V583_genome.tsv
#blastn -query ${file} -db ${PLASMID_FINDER}/plasmidfinder_db/AllRep -perc_identity 90 -num_threads 16 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/CoAsm_rep.tsv
#blastn -query ${file} -db ${HIP}/HIP11704 -perc_identity 90 -num_threads 16 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/CoAsm_HIP_genome.tsv
blastn -query ${file} -db ${PLS} -perc_identity 90 -num_threads 16 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/CoAsm_plsdb.tsv



#for f in ${ASM}/blast_result/*.tsv;do
#	SAMPLE=${file##*/}
#	BASE=${SAMPLE%.tsv}
#	awk '$8 >=800' ${file} > ${ASM}/blast_result/filtered/${BASE}_800_90.tsv
#done


#>NC_004668.1 Enterococcus faecalis V583 chromosome, complete genome
#>NC_004669.1 Enterococcus faecalis V583 plasmid pTEF1, complete sequence
#>NC_004671.1 Enterococcus faecalis V583 plasmid pTEF2, complete sequence
#>NC_004670.1 Enterococcus faecalis V583 plasmid pTEF3, complete sequence


