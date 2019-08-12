####New 150PE reads

1. individual assembly
2. Co-assembly: bin3C, combine HiC reads and get bam before subjected to bin3C
3. CRISPR array in co-assembly

####Co-assembly binning 

Script: /net/ab/cb/68/cxz163430/project/HiC/script/coassembly_150pe.sh

1. update metwrap on aleph & fix the endless warning with concoct
2. metawrap binning with metaBat2, Maxbin and CONCOCT

```sh
gunzip ${prc_data}/SG_Reads1_combined_150.fq.gz
gunzip ${prc_data}/SG_Reads2_combined_150.fq.gz

mv ${prc_data}/SG_Reads1_combined_150.fq ${prc_data}/SG_combined_150_1.fastq
mv ${prc_data}/SG_Reads2_combined_150.fq ${prc_data}/SG_combined_150_2.fastq

source ~/.bashrc
unset PYTHONPATH
conda activate metawrap-env
metawrap binning -o ${BIN}/Metawrap -t 24 -a ${ASSM}/final.contigs.fa --metabat2 --maxbin2 --concoct ${prc_data}/SG_combined*fastq

conda deactivate
```
3. bin3c

```sh
bwa index ${ASSM}/final.contigs.fa
bwa mem -t 24 -5SP ${ASSM}/final.contigs.fa  '<zcat ${HCdat}/*1_paied.fastq.gz' '<zcat  ${HCdat}/*2_paired.fastq.gz'  | \
    samtools view -h -F 2316 -bS -@ 24 - | \
    samtools sort -n -o ${HCdat}/Coassem_150pe_2316.bam -
##download to local work station, cannot let bin3C work on either mz or aleph....

./bin3C.py mkmap -e MluCI -e Sau3AI -v /ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/final.contigs.fa /ms/11/cong/project/HiC/BINNING/150PE_coassembly/ /ms/11/cong/project/HiC/BINNING/150PE_coassembly/bin3C

```

####Individual assembly

Script:

1. bin3C binning on local work station (150pe_bin3c.sh)

	results stored at /ms/11/cong/project/HiC/BINNING/150PE/bin3c

	upload bins with size > 15000bp to mz


2. 