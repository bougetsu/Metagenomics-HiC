#### Things to do

1. check Vancomycin location of pTEF plasmid and the changes caused by plasmid introduction.

    - pTEF1 and pTEF2 are structurally similar to the archetypal pheromone responsive plasmids pAD1 (15) and pCF10 (14), respectively, and pTEF3 belongs to the family of pAMβ1 broad host range plasmids.
    - The unique regions in pTEF1 include a Tn4001-like transposon encoding aminoglycoside resistance and another IS-flanked element carrying erythromycin resistance and multidrug resistance genes.
    - V583: VanB-Tn1549

2. Complare the plasmid linked E.f contigs in 3 samples (it is wired that VE has more plasmid links)

    - more plasmid contigs in VE???
    - plasmid contig abundance in three samples (TPM, normalized by total reads?)
    - Hi-C links of contigs - bin in three samples
    - Bin abundance


3. Check the insertion genes (flanking genes) of resistance genes, if the sequence changed over time?? usingi individual assembly


ClobberError: This transaction has incompatible packages due to a shared path.
  packages: defaults::matplotlib-2.2.3-py27hb69df0a_0, conda-forge::matplotlib-base-2.2.4-py27hfd891ef_0
  path: 'lib/python2.7/site-packages/pylab.pyc'

4. prepare replicon and resistance gene region on contigs in bed format (coassembly)



`samtools view -b -L region.bed test.bam > test_region.bam`

sort and index each bam file, and then extract rep region

####Day 13, HIP

| contig       | Links_number | ref           | replicon |
|--------------|--------------|---------------|----------|
| k141_813771  | 1            | pTEF1+VE14089 | Rep9     |
| k141_2242115 | 13           | pTEF1         | Rep1     |
| k141_237538  | 471          | HIP_pls       | Rep1     |
| k141_699651  | 149          | gb|GG692639.1 | Rep9     |




check the reads mapped on to pTEF1, maybe can be mapped on to HIP replicon

1. reads name, sequence mapped on to HIP replicon?? why bwa assign then to pTEF1??

blast??

```
HIP='/ms/11/cong/project/HiC/data/ref/HIP11704'
ASM='/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE'
blastn -query ${file} -db ${HIP}/HIP11704 -perc_identity 90 -num_threads 16 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/CoAsm_HIP_genome.tsv

PLASMID_FINDER='/home/cxz163430/software/plasmidfinder'
blastn -query ${file} -db ${PLASMID_FINDER}/plasmidfinder_db/AllRep -perc_identity 90 -num_threads 16 -outfmt '6 qseqid sseqid pident mismatch gapopen qlen slen length qstart qend sstart send evalue' -out ${ASM}/blast_result/CoAsm_rep.tsv
```


* Replicon: No hit
* HIP genome: K00364:132:H2VWMBBXY:4:1120:2696:40631    gb|GG692631.1|, rep9 region?


* very tedious thing!!!

match seq on reads K00364:132:H2VWMBBXY:4:2202:32086:2615@k141_2242115: GATCGATCATTAAAGTAAAAAAAGAAGAAAAAGAAAGCTATATAAAGG

can be also mapped to HIP_PLS contigs (k141_237538)



So we only keep those reads on HIP_pls and HIP_chr replicon region.

####Day 13, V583

####Day 13, VE
Rep9 on VE chr :3094878 to 3095888
rep9 on pTEF2 ref seq:100 to 1100
Clostridioides difficile strain DSM 27639
Clostridium difficile transposon Tn6218

k141_237538


2. contig abundance


Find res gene on individual assembly and the prokka annotate those res contigs.

* resfinder on mz
* get res contigs in R
*
seqkit grep -f resfinder/0HIP11704_contigs.list 0HIP11704_150PE.contigs.fa > resfinder/0HIP11704_res.contigs.fa
prokka --quiet --cpus $threads --outdir ${out}/prokka_out/$bin_name --prefix $bin_name ${out}/tmp_bin.fa
* PROKKA




#####VE STRAIN

2011 PAPER

"Large-Scale Screening of a Targeted Enterococcus faecalis Mutant Library Identifies Envelope Fitness Factors"

This strain has not undergone major DNA rearrangements, as indicated by DNA microarray analysis, pulsed field gel electrophoresis, Southern blotting, long-range PCR and DNA sequencing analyses, with one exception: a 20.5-kb region of pTEF1 (from efa0063 up to efa0006) is integrated between chromosomal genes ef3209 and ef3210, and is flanked by two copies of IS1216.

2019 PAPER
we used a plasmid-cured
derivative of strain V583 named VE14089 that was obtained after novobiocin and
thermic shocks (22). Strain VE14089 was generated from our clone of the V583 strain
that was renamed VE14002 and that lacks the pTEF3 plasmid (23).




```
python3 ~/software/SolidBin/SolidBin.py --contig_file assembly.fa  --coverage_profiles concoct_depth_local.txt --composition_profiles kmer_4_tmp.csv  --output test_solidbin/result.csv
```


####did HIP pls and VanHAX shared same transfering pattern?


1. restricted onto replicon and res gene region
  extracted linked contig and get their belonged bins

  a. use rep and res region: only bin.19 and bin.141 while a few NA


```{sh}
samtools view -b 13HIP11704_hic2CoAsm2316_150pe.sorted.bam k141_237538 > 13HIP11704_hic2CoAsm2316_150pe.k141_237538.bam
```


2. taxonomic classification of all linked contigs

  a. get contig names from links file and then get fasta.

  b. using blobtools

  ```{sh}
  conda activate blobtools
  export ASM='/net/ab/cb/68/cxz163430/project/HiC/ASSEMBLY/Megahit_coasm_150pe'
  MAP='/net/ab/cb/68/cxz163430/project/HiC/data/mapping_150CoAsm'
  ./blobtools create \
 -i ${ASM}/final.contigs.fa \
 -b ${MAP}/0HIP11704_150PE.bam \
 -b ${MAP}/0V583_150PE.bam \
 -b ${MAP}/0VE14089_150PE.bam \
 -b ${MAP}/13HIP11704_150PE.bam \
 -b ${MAP}/13V583_150PE.bam \
 -b ${MAP}/13VE14089_150PE.bam \
 -t ${ASM}/assembly.vs.nt.mts1.hsp1.1e25.megablast.out \
 -o ${ASM}/blobplot_nt

  ```


ERROR:
Traceback (most recent call last):
  File "/net/ab/cb/68/cxz163430/software/blobtools/blobtools", line 7, in <module>
    main()
  File "/net/ab/cb/68/cxz163430/software/blobtools/lib/interface.py", line 60, in main
    create.main()
  File "/net/ab/cb/68/cxz163430/software/blobtools/lib/create.py", line 119, in main
    blobDb.parseCoverage(covLibObjs=cov_libs, estimate_cov=estimate_cov_flag, prefix=prefix)
  File "/net/ab/cb/68/cxz163430/software/blobtools/lib/BtCore.py", line 369, in parseCoverage
    base_cov_dict, covLib.reads_total, covLib.reads_mapped, read_cov_dict = BtIO.parseBam(covLib.f, set(self.dict_of_blobs), estimate_cov)
  File "/net/ab/cb/68/cxz163430/software/blobtools/lib/BtIO.py", line 224, in parseBam
    base_cov_dict, read_cov_dict = estimate_coverage(aln, set_of_blobs)
  File "/net/ab/cb/68/cxz163430/software/blobtools/lib/BtIO.py", line 232, in estimate_coverage
    est_read_length = estimate_read_lengths(aln, set_of_blobs)
  File "/net/ab/cb/68/cxz163430/software/blobtools/lib/BtIO.py", line 201, in estimate_read_lengths
    for read in aln.fetch(header):
  File "pysam/libcalignmentfile.pyx", line 1081, in pysam.libcalignmentfile.AlignmentFile.fetch
  File "pysam/libchtslib.pyx", line 686, in pysam.libchtslib.HTSFile.parse_region
ValueError: invalid contig `k141_1531323`



```{sh}
diamond blastx \
 --query final.contigs.fa \
 --db ~/data/NCBI_nr/nr.dmnd \
 --outfmt 6 \
 --sensitive \
 --max-target-seqs 1 \
 --evalue 1e-25 \
 --out assembly.vs.nr.1e25.diamond.out


diamond blastx  --threads 36 --query final.contigs.fa  --db ~/data/NCBI_nr/nr.dmnd --outfmt 6  --sensitive  --max-target-seqs 1  --evalue 1e-25  --out assembly.vs.nr.1e25.diamond.out

```



```{sh}
export ASM='/net/ab/cb/68/cxz163430/project/HiC/ASSEMBLY/Megahit_coasm_150pe'
krakenuniq --db ~/data/Krakenuniq_db --threads 24 ${ASM}/final.contigs.fa --report-file ${ASM}/final.contigs.krareport > ${ASM}/final.contigs.kra

```



3. recheck the heatmap of 13HIP, the links from HIP_pls. The # of links in heatmap is (2*links from bam file)/bin abundance.

4. Check the taxonomic label of all the contigs linked with HIP_pls in 13HIP

"~/Documents/Work/project/HiC/Result/13HIP_pls_linked.kraken"




5. For those "E.f" contigs in other bins, manually remove them and add back to E.f bin, check the completedness and contamination rate



6. Use salmon to calculate the  abundance of all contigs across the samples
