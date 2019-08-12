#! /bin/bash
echo TRIM START

RAW_READS='/ms/11/cong/data/raw_data/HiC'
prc_data='/ms/11/cong/project/HiC/processed_reads'
TRIMM='/home/cxz163430/software/Trimmomatic-0.38/trimmomatic-0.38.jar'
ASSM='/ms/11/cong/project/HiC/ASSEMBLY'

mkdir -p ${prc_data}
mkdir -p ${ASSM}/IDBAUD
mkdir ${ASSM}/MEGAHIT
for F in ${RAW_READS}/*_R1_001.fastq.gz; do 
	R=${F%R1*}R2_001.fastq.gz
	BASE=${F##*/}
	SAMPLE=${BASE%_L00*}
	echo $F
	echo $R
	echo $SAMPLE
	#bbduk.sh -Xmx80g  in1=$F in2=$R ref=/home/cxz163430/software/bbmap/resources/adapters.fa out1=${prc_data}/${SAMPLE}_1.fastq.gz out2=${prc_data}/${SAMPLE}_2.fastq.gz k=23 ktrim=r mink=11 hdist=1 minlength=50 tpe tbo ftm=5 qtrim=r trimq=10

	if [[ "$SAMPLE" =~ "_SG_" ]];then
		echo "##################################"
		echo "MEGAHIT"
		echo $SAMPLE
		echo "##################################"
		megahit -1 ${prc_data}/${SAMPLE}_1.fastq.gz -2 ${prc_data}/${SAMPLE}_2.fastq.gz -o ${ASSM}/MEGAHIT/${SAMPLE}
	fi
done	

