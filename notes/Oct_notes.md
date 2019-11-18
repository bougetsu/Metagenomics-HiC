
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


  ```

  ERROR: invalid contig "k141_605712".

  found in header of 13 days bam files but not day 0.

  solution to try:

  1. blobtools within metawrap, use bowtie2 to get bam file and then re-do

    1.a. try 0HIP.bowtie.bam on ab1

    1.b. try 0HIP.bam on ab1

  2. individual blobtools with 13HIP (ab1)...not working, "ValueError: invalid contig `k141_841324`"
  3. try just taxify



  4. check the rep_res contigs, alignment info (region)

  5. plot the linked contigs by samples




3. recheck the heatmap of 13HIP, the links from HIP_pls. The # of links in heatmap is (2*links from bam file)/bin abundance.

4. Check the taxonomic label of all the contigs linked with HIP_pls in 13HIP

"~/Documents/Work/project/HiC/Result/13HIP_pls_linked.kraken"





5. For those "E.f" contigs in other bins, manually remove them and add back to E.f bin, check the completedness and contamination rate

  a. criteria: classified as "E.f" by kraken and refpro but not in the bin.19
  b. get all contigs labeled by Kraken, ref_pro and nt, get those not in bin.19 but "grep with Enterococcus"
  c. in nt, some contigs may be assgined with different labels and similar scores
  `"ASSEMBLY/CoAsm_150PE/Taxify/metawrap/nt_confused_contig.list"`
  c. give a flag to it: Y : both E.f in kraken and nt, N: not E.f in nt and ref_pro, F: E.fm in kraken and nt

    result in `ASSEMBLY/CoAsm_150PE/Taxify/metawrap/To_check_ef.tsv`

    only several binned contigs was labeled as "Y", since most of them have in consistent labels for kra, nt and ref_pro. How to do?? some need to move out from original bins.





  d. got those has flag == Y and add to bin.19, use checkM to evaluate the completeness and contamination rate

    `"/ms/11/cong/project/HiC/ASSEMBLY/CoAsm_150PE/Taxify/metawrap/new_ef_contig_v1.list"`



6. Use salmon to calculate the  abundance of all contigs across the samples
