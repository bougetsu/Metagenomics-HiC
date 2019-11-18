####


#### steps

* mapping Hi-C reads from each sample to co-assembly

* Plasmid contigs:
    1. blast plasmid against V583, HIP genome, replicon and PLS-DB (plasmid database)

        1. filter criteria: identity > 90, length > 800 bp, matched length/ contig length > 0.6
    2. candidate contig results: blast against nr data base using webpage blast and manually checked the results, removed those having hits on other genomes
    3. Results: 29 contigs with 6 carrying Rep
* Res contig:
    1. using ResFinder

* check HiC links in each sample by plasmid or ARG       
    1. find some contigs connected with EF but not clustered in dRep results, gp back to bin3C results, many of them are in EF cluster by bin3C results
    2. add thost contigs back into EF bin (linked with ef plasmid > 4 and binned into EF by bin3c)
    3. ~ 80 bins in EF in dRep and 208 in EF in bin3c, check completedness and contamination rate before and after adding contigs from bin3C

| bin    | size | completedness | contamination rate | Strain heterogeneity |
|--------|------|---------------|--------------------|----------------------|
| before | 82   | 98.41         | 0.37               | 50                   |
| after  | 141  | 99.44         | 0.94               | 50.00                |

* HiC links summarized


|                | V583    |         | HIP     |         | VE      |         |
|----------------|---------|---------|---------|---------|---------|---------|
| Day            | 0       | 13      | 0       | 13      | 0       | 13      |
| total links    | 2998446 | 5364364 | 2655326 | 9029966 | 2612002 | 4363302 |
| non-self links | 2789944 | 4513821 | 2447986 | 8371705 | 2430636 | 3690062 |
| plasmid links  | 4       | 174     | 0       | 1178    | 2       | 1688    |
| non clustered  | 2       | 32      | 0       | 362     | 2       | 318     |
| ARG links      | 484     | 936     | 738     | 5352    | 468     | 986     |
| non clustered  | 302     | 356     | 374     | 2990    | 232     | 374     |


* draw heatmap showing HiC links between bins and plasmid/ARG contigs
    1. Hi-C links were normalized by the (absolute) abundance of the bin (genome) Transcripts Per Million ? or NumReads?
    2.  how about length?? should it be length * abundance = reads??




#### To check:

1. VanHBX on V583 chr

2. VanHAX on HIP Plasmid

3. - ANT(6)-Ib is an aminoglycoside nucleotidyltransferase gene encoded by transferable pathogenicity islands in C. fetus subsp. fetus and B. subtilis

  - Clostridium perfringen

4. ErmB in V583: Plasmid pTEF1
k141_1233618

5. blast the resistance contigs and find their origins/host

6. Flanking region of resistance genes:

  - PROKKA, blast, compare with plasmid, sequence divergence, align those genes

      * using metawrap ANNOTATE_BIN on alphe ab1

      * got ptoblem with perl5 lib, commented those lines which manually set perlib path, worked!

      * for those un-clustered contig... try prokka with all contigs?
  - use individual assembly to check
    * PROKKA on mz ... = = canot let it work...
  - how to compare????  hard to think

7. Got VE14089 strain from Sara on Sep 11th.

  - streak BHI plate on Sep 11th
  - pick colonies (not mono-colony in case of mixed/contaminated sample) and culture overnight to make 3 stocks
  - sequence, talk with luke about that. extract DNA and get PO
    * will get email and protocal from dennise
  - PCR (later)

  - streak on 5 mg/ml Van to see if it can grow.
