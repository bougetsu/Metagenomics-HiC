####New 150PE reads

1. individual assembly
2. Co-assembly: bin3C, combine HiC reads and get bam before subjected to bin3C
3. CRISPR array in co-assembly

####Co-assembly binning 

Script: /net/ab/cb/68/cxz163430/project/HiC/script/coassembly_150pe.sh

1. update metwrap on aleph & fix the endless warning with concoct
2. metawrap binning with metaBat2, Maxbin and CONCOCT

```bash
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

get error msg

```sh
------------------------------------------------------------------------------------------------------------------------
-----             looks like our default perl libraries are not the conda ones. Manually setting perl5             -----
-----                                              library directory                                               -----
------------------------------------------------------------------------------------------------------------------------

metawrap path: /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/bin/metawrap
Will use perl5 libraries located in /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/lib/perl5/site_perl/5.22.0 - hopefully they are there...

------------------------------------------------------------------------------------------------------------------------
-----                                       Starting binning with MaxBin2...                                       -----
------------------------------------------------------------------------------------------------------------------------

Can't locate HTTP/Status.pm in @INC (you may need to install the HTTP::Status module) (@INC contains: /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/lib/perl5/site_perl/5.22.0/x86_64-linux-thread-multi /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/lib/perl5/site_perl/5.22.0 /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/lib/site_perl/5.26.2/x86_64-linux-thread-multi /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/lib/site_perl/5.26.2 /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/lib/5.26.2/x86_64-linux-thread-multi /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/lib/5.26.2 .) at /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/lib/site_perl/5.26.2/LWP/Simple.pm line 14.
BEGIN failed--compilation aborted at /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/lib/site_perl/5.26.2/LWP/Simple.pm line 14.
Compilation failed in require at /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/bin/run_MaxBin.pl line 4.
BEGIN failed--compilation aborted at /net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/bin/run_MaxBin.pl line 4.

```


To install that
```
conda install -c bioconda perl-http-message 
```


3. bin3c

```bash
bwa index ${ASSM}/final.contigs.fa
bwa mem -t 24 -5SP ${ASSM}/final.contigs.fa  '<zcat ${HCdat}/*1_paied.fastq.gz' '<zcat  ${HCdat}/*2_paired.fastq.gz'  | \
    samtools view -h -F 2316 -bS -@ 24 - | \
    samtools sort -n -o ${HCdat}/Coassem_150pe_2316.bam -
##download to local work station, cannot let bin3C work on either mz or aleph....

./bin3C.py mkmap -e MluCI -e Sau3AI  --min-reflen 1500 -v /ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/final.contigs.fa /ms/11/cong/project/HiC/BINNING/150PE_coassembly/Coassem_150pe_2316.bam /ms/11/cong/project/HiC/BINNING/150PE_coassembly/bin3C_1k5

bin3C  cluster --no-spades --only-large  -v /ms/11/cong/project/HiC/BINNING/150PE_coassembly/bin3C/contact_map.p.gz /home/cxz163430/Project/HiC/bin3C/cluster


###USING DE NISCO
bin3C  cluster --no-spades --only-large  -v /home/cxz163430/Project/HiC/data/contact_map.p.gz /home/cxz163430/Project/HiC/bin3C/cluster

bin/python2 ./bin3C.py mkmap -e MluCI -e Sau3AI  --min-reflen 1500 -v /home/cxz163430/Project/HiC/data/CoAssembly_150pe.contigs.fa /home/cxz163430/Project/HiC/data/Coassem_150pe_2316.bam /home/cxz163430/Project/HiC/bin3C/map_1k5



```

__Got memory error__ : bin3C failed to run on mz and aleph, and my desktop doesn't have enough memmory :)


3. blobtools

```bash
blastn \
 -query ~/project/HiC/ASSEMBLY/Megahit_coasm_150pe/final.contigs.fa \
 -db ~/data/NCBI_nt/nt \
 -outfmt ’6 qseqid staxids bitscore std’ \
 -max_target_seqs 10 \
 -max_hsps 1 \
 -evalue 1e-25 \
 -out ~/project/HiC/ASSEMBLY/Megahit_coasm_150pe/assembly.vs.nt.mts1.hsp1.1e25.megablast.out

diamond blastx \
 --query $ASSEMBLY \
 --db uniprot_ref_proteomes.diamond.dmnd \
 --outfmt 6 \
 --sensitive \
 --max-target-seqs 1 \
 --evalue 1e-25 




blobtools taxify \ 
 -f diamond.out \
 -m uniprot_ref_proteomes.taxids 
 -s 0 \ # column of sequenceID of subject in taxID mapping file
 -t 2 # column of TaxID of sequenceID in taxID mapping file
```



####Individual assembly

Script:

1. bin3C binning on local work station (150pe_bin3c.sh)

	results stored at /ms/11/cong/project/HiC/BINNING/150PE/bin3c

	upload bins with size > 15000bp to mz


2. metawrap bining

got error with mz36, 38, 40: 

```sh
------------------------------------------------------------------------------------------------------------------------
-----                                  Running CheckM on maxbin2_bins.contigs bins                                 -----
------------------------------------------------------------------------------------------------------------------------


*******************************************************************************
 [CheckM - tree] Placing bins in reference genome tree.
*******************************************************************************


Unexpected error: <type 'exceptions.OSError'>
Traceback (most recent call last):
  File "/net/zmf2/cb/9/cxz163430/miniconda3/envs/metawrap-env/bin/checkm", line 708, in <module>
    checkmParser.parseOptions(args)
  File "/net/zmf2/cb/9/cxz163430/miniconda3/envs/metawrap-env/lib/python2.7/site-packages/checkm/main.py", line 1251, in parseOptions
    self.tree(options)
  File "/net/zmf2/cb/9/cxz163430/miniconda3/envs/metawrap-env/lib/python2.7/site-packages/checkm/main.py", line 110, in tree
    binFiles = self.binFiles(options.bin_folder, options.extension)
  File "/net/zmf2/cb/9/cxz163430/miniconda3/envs/metawrap-env/lib/python2.7/site-packages/checkm/main.py", line 87, in binFiles
    all_files = os.listdir(binFolder)
OSError: [Errno 20] Not a directory: 'maxbin2_bins.contigs'

```

Try to re-run:

mz 38 and 40 finished binning successfully, while 36 still got the same error :)




####Anvio

##### set up server

got problems.....= =
need to re-install




###Reads levels

got 2308 and 2312 HiC reads

1. Kraken

There are about 1/4 ~ 1/3 of 2312 reads cannot be classified by kraken

2. Diamond

Is it right to using nr?? what if the reads mapped to a un-translated region??


###CRISPR




###REMOVE PYTHON AND CONDA on alepha

###re-install python2
http://thelazylog.com/install-python-as-local-user-on-linux/

```sh
mkdir /net/ab/cb/68/cxz163430/local/python2.7
cd python2.7
wget "https://www.python.org/ftp/python/2.7.16/Python-2.7.16.tgz"
tar xzvf Python-2.7.16.tgz
find /net/ab/cb/68/cxz163430/local/python2.7 -type d | xargs chmod 0755
cd Python2.7.16
./configure --prefix=/net/ab/cb/68/cxz163430/local/python2.7
make && make install
cd ..
wget https://bootstrap.pypa.io/get-pip.py && python get-pip.py --user


```

```Executing transaction: - WARNING conda.core.envs_manager:register_env(46): Unable to register environment. Path not writable or missing.

  environment location: /net/ab/cb/68/cxz163430/miniconda3

  registry file: /net/ab/cb/68/cxz163430/.conda/environments.txt

done

installation finished.

WARNING:
    You currently have a PYTHONPATH environment variable set. This may cause
    unexpected behavior when running the Python interpreter in Miniconda3.
    For best results, please verify that your PYTHONPATH only points to
    directories of packages that are compatible with the Python interpreter
    in Miniconda3: /net/ab/cb/68/cxz163430/miniconda3
```



###Individual mapping bam files to co-assembly

```
##running on mz
###bwa

ASM='/net/zmf2/cb/9/cxz163430/project/Hic/Assembly/CoAssembly_150PE'
HiCMap='/net/zmf2/cb/9/cxz163430/project/Hic/Binning_150PE/CoAsm/HiC_mapping'
prc_data='/net/zmf2/cb/9/cxz163430/project/Hic/processed_data'

bwa mem -t 8 -5SP ${ASM}/final.contigs.fa ${prc_data}/${SAMPLE}_1_paied.fastq.gz ${prc_data}/${SAMPLE}_2_paired.fastq.gz | \
    samtools view -h -F 2316 -bS -@ 24 - | \
    samtools sort -n -o ${HiCMap}/${BASE}_hic2CoAsm2316_150pe.bam -

##check total HiC counts

##Normalize it?? across samples


```


/net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/bin/metawrap-scripts/split_salmon_out_into_bins.py /net/ab/cb/68/cxz163430/project/HiC/Binning/BINNING_CO_150PE/QUANT_BINS_RF/quant_files/ /net/ab/cb/68/cxz163430/project/HiC/Binning/BINNING_CO_150PE/BIN_REFINEMENT/metawrap_50_10_bins /net/ab/cb/68/cxz163430/project/HiC/ASSEMBLY/Megahit_coasm_150pe/final.contigs.fa  > /net/ab/cb/68/cxz163430/project/HiC/Binning/BINNING_CO_150PE/QUANT_BINS_RF/bin_abundance_table.tab


/net/ab/cb/68/cxz163430/miniconda3/envs/metawrap-env/bin/metawrap-scripts/make_heatmap.py /net/ab/cb/68/cxz163430/project/HiC/Binning/BINNING_CO_150PE/QUANT_BINS_RF/bin_abundance_table.tab /net/ab/cb/68/cxz163430/project/HiC/Binning/BINNING_CO_150PE/QUANT_BINS_RF/bin_abundance_heatmap.png

##de-replicate binning results.. try it later

need seperate bining results and coassembly


```
~/software/drep/bin/dRep dereplicate metawrap -g metawrap_50_10_bins/*.fa 
```


##get CRISPR and Res/Rep results


```
/home/cxz163430/software/CRISPRCasFinder//home/cxz163430/software/CRISPRCasFinder

perl  -keep -i /ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/final.contigs.fa -o /ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/crisper -html -gscf -cas -cpuM 16
```


ARG homologues were identified using BLASTN with the nucleotide sequences extracted from the Prodigal ORF locations as a query against the transferrable ARG ResFinder database [57]. Hits with a minimum 95% nucleotide sequence identity and 90% ARG sequence coverage were retained as candidate ARGs. 



```

conda activate blobtools
export PYTHONPATH=/net/ab/cb/68/cxz163430/miniconda3/lib/python3.7/site-packages
bblobtools taxify 
ASM='/net/ab/cb/68/cxz163430/project/HiC/ASSEMBLY/Megahit_coasm_150pe'
blobtools taxify -f ${ASM}/assembly.vs.refpro.1e25.diamond.out -m /net/ab/cb/68/cxz163430/data/Uniprot/uniprot_ref_proteomes.taxids -s 0 -t 2 -o ${ASM}/blob_taxify



###ab1
python3 resfinder.py -i ${ASM}/final.contigs.fa -o ${ASM}/ResFinder -p /net/ab/cb/68/cxz163430/software/resfinder_db -mp ~/local/blast/bin/blastn -t 0.90 -l 0.60 > ${ASM}/resfinder.log


##local CRIPRFINDER

##ab2

```

1. blobtools

2. plots showing the abundance of each sample using taxo annotated bins

3. summary of each HiC samples aligned to coassembly



```
dat='/ms/11/cong/data/raw_data/HiC'
out='/ms/11/cong/project/HiC/data/ReadsQC/HiC'
for F in ${dat}/*HC*R1_001.fastq.gz; do 
  R=${F%R1_001.fastq.gz}R2_001.fastq.gz
  BASE=${F##*/}
  SAMPLE=${BASE%PALMER*}
  fastqc -q -t 12 -o ${out} -f fastq $F $R
done  

fastqc -q -t $threads -o ${out}/pre-QC_report -f fastq $reads_1 $reads_2

```




4. get HiC coverage of each samples on coassembly
5. Prodigal

```
dat='/net/ab/cb/68/cxz163430/project/HiC/Binning/BINNING_CO_150PE/BIN_REFINEMENT'
for F in ${dat}/metawrap_50_10_bins/*.fa; do 
  BASE=${F##*/} 
  id=${BASE%.fa} 
  prodigal -p meta -a ${dat}/proteins/${id}.faa -q -i ${F} -f gff -o ${dat}/proteins/${id}_prodigal.gff  
done  

```

SOurmash

`/net/ab/cb/68/cxz163430/data/sourmash/genbank-d2-k31.sbt.json`

```
sourmashgb='/ms/11/cong/data/Sourmash/genbank-d2-k31.sbt.json'
dat='/ms/11/cong/project/HiC/BINNING/150PE_coassembly/BIN_REFINEMENT/metawrap_dRep'
for F in ${dat}/dereplicated_genomes/*.fa; do 
  BASE=${F##*/} 
  id=${BASE%.fa}
  sourmash compute --scaled 1000 -k 31 -o ${dat}/sourmash/${id}.sig $F 
  sourmash gather -k 31 ${dat}/sourmash/${id}.sig ${sourmashgb} -o ${dat}/sourmash/${id}.csv > ${dat}/sourmash/${id}.sm
done 

```


6. dRep

202 genomes were reduced to 162.

Do annotation and quant_bin at aleph

```
/net/ab/cb/68/cxz163430/project/HiC/Binning/BINNING_CO_150PE/metawrap_dRep/dereplicated_genomes
```


7. Get Hi-C links using connects.R

  a. contigs class: plasmid, (rep), ARG

  confirm plasmid contigs and ARG

  using BLASTn results (V583, HIP, rep, PLS-DB) get plasmid contig if list, extract contig fasta to blast using website to confirm whether unique to plasmid sequence all also found hits on chr

  ```
  #/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE
  seqkit  grep -f blast_result/plasdmid_contig.list final.contigs.fa > blast_result/plasdmid_contig.fa
  ```

  save manually checked results to `plasmid_contig_checked.info`

######
check k141_607023 cluster in bin3c and which is ef cluster in bin3c

~ 80 bins in bin.19 in dRep and 208 in bin.38 in bin3c

check completedness and contamination rate in bin.19+ other contigs

after adding contigs from bin3C

| bin    | size | completedness | contamination rate | Strain heterogeneity |
|--------|------|---------------|--------------------|----------------------|
| before | 82   | 98.41         | 0.37               | 50                   |
| after  | 141  | 99.44         | 0.94               | 50.00                |

Links STAT

|                | V583    |         | HIP     |         | VE      |         |
|----------------|---------|---------|---------|---------|---------|---------|
| Day            | 0       | 13      | 0       | 13      | 0       | 13      |
| total links    | 2998446 | 5364364 | 2655326 | 9029966 | 2612002 | 4363302 |
| non-self links | 2789944 | 4513821 | 2447986 | 8371705 | 2430636 | 3690062 |
| plasmid links  | 4       | 174     | 0       | 1178    | 2       | 1688    |
| non clustered  | 2       | 32      | 0       | 362     | 2       | 318     |
| ARG links      | 484     | 936     | 738     | 5352    | 468     | 986     |
| non clustered  | 302     | 356     | 374     | 2990    | 232     | 374     |

