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

./bin3C.py mkmap -e MluCI -e Sau3AI -v /ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/final.contigs.fa /ms/11/cong/project/HiC/BINNING/150PE_coassembly/Coassem_150pe_2316.bam /ms/11/cong/project/HiC/BINNING/150PE_coassembly/bin3C

bin3C  cluster --no-spades --only-large  -v /ms/11/cong/project/HiC/BINNING/150PE_coassembly/bin3C/contact_map.p.gz /home/cxz163430/Project/HiC/bin3C/cluster


###USING DE NISCO
bin3C  cluster --no-spades --only-large  -v /home/cxz163430/Project/HiC/data/contact_map.p.gz /ms/11/cong/project/HiC/BINNING/150PE_coassembly/bin3C/cluster
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