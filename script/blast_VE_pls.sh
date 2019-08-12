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

ASM='/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT'
RES=''
for file in ${ASM}/*VE*; do
	BASE=${file##*/}
	echo $BASE
	echo ${file}
	blastn -query ${file}/final.contigs.fa -db ${V583}/V583 -perc_identity 90 -num_threads 8  -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/${BASE}_genome.tsv
	blastn -query ${file}/final.contigs.fa -db ${PLASMID_FINDER}/plasmidfinder_db/AllRep -perc_identity 90 -num_threads 8 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/${BASE}_rep.tsv
	blastn -query ${file}/final.contigs.fa -db ${PLS} -perc_identity 90 -num_threads 8 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/${BASE}_plsdb.tsv
done	

done

###replicon
