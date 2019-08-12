
#! /bin/bash
echo START

Dat='/net/ab/cb/68/cxz163430/project/HiC/data/150PE'
prc_data='/net/ab/cb/68/cxz163430/project/HiC/data'
ASSM='/net/ab/cb/68/cxz163430/project/HiC/ASSEMBLY/Megahit_coasm_150pe'
BIN='/net/ab/cb/68/cxz163430/project/HiC/'


for f in ${Dat}/*_R1*;do
	R=${F%_R1*}_R2_001.fastq.gz
	SAMPLE=${F##*/}
	BASE=${SAMPLE%%-*}
	echo ${BASE}
	bbduk.sh -Xmx60g  in1=${F} in2=${R} ref=/net/zmf2/cb/9/cxz163430/software/bbmap/resources/adapters.fa out1=${prc_data}/${BASE}_150PE_R1.fastq.gz out2=${prc_data}/${BASE}_150PE_R2.fastq.gz k=23 ktrim=r mink=11 hdist=1 minlength=50 tpe tbo ftm=5 qtrim=r trimq=10
done



#mkdir -p ${ASSM}
cd ${prc_data}
cat *_1.fastq.gz > SG_Reads1_combined.fq.gz
cat *_2.fastq.gz > SG_Reads2_combined.fq.gz
#echo START_MEGAHIT
megahit --kmin-1pass -1 ${prc_data}/SG_Reads1_combined.fq.gz -2 ${prc_data}/SG_Reads2_combined.fq.gz -o ${ASSM}/megahit_coasm

#cd ${BIN}/bin3c/clust_default/fasta/
#rename CL bin. ${BIN}/bin3c/clust_default/fasta/*
#rename fna fa ${BIN}/bin3c/clust_default/fasta/*
#sed -i "s#>CL[0-9]*_[0-9]* contig:\(k[0-9]*_[0-9]*\) ori:UNKNOWN length:[0-9]*#>\1#g" ${BIN}/bin3c/clust_default/fasta/*
#cd $HOME/project/HiC/script

#echo #####RENAME FINISHED#######

#source ~/.bashrc
#conda activate metawrap-env
#metawrap binning -o INITIAL_BINNING_CO_0423 -t 12 -a  ${ASSM}/megahit_coasm/final.contigs.fa --metabat2 --maxbin2 --concoct ${prc_data}/*fastq
#mkdir -p ${BIN}/REFINEMENT_5k
#metawrap bin_refinement -o ${BIN}/REFINEMENT_5k -t 10 -A ${BIN}/metabat2_bins/ -B ${BIN}/maxbin2_bins/ -C ${BIN}/bin3c/clust_default/5k_fasta/  -c 50 -x 10



#gunzip -c ${prc_data}/PooledReads/SG_Reads1_combined.fq.gz > ${prc_data}/PooledReads/SG_combined_1.fastq
#gunzip -c ${prc_data}/PooledReads/SG_Reads2_combined.fq.gz > ${prc_data}/PooledReads/SG_combined_2.fastq

#echo "checkM on BLOB"
#source ~/.bashrc
#unset PYTHONPATH
#conda activate metawrap-env
#metawrap blobology -a ${ASSM}/megahit_coasm/final.contigs.fa -t 4 -o ${BIN}/BLOBOLOGY_mw --bins ${BIN}/REFINEMENT_5k/metawrap_50_10_bins ${prc_data}/PooledReads/SG_combined_*.fastq
#metawrap quant_bins -b ${BIN}/REFINEMENT_5k/metawrap_50_10_bins -o ${BIN}/QUANT_BINS -a ${ASSM}/megahit_coasm/final.contigs.fa ${prc_data}/
#echo "##################################"
#echo "##################################"

#conda deactivate

