
####reformat HiC counts as Qiwei suggested

Row: contigs:plasimd (pre-/post blast filter) & resistance
Col: linked contigs with taxonomic label (kraken)

a. check the tax label between kraken and metawrap for contigs
  how to deal with inconsistency

    use kraken for now

b. take the abundance into consideration??

    since qiwei needs counts, provide depth file to him separately

2,387,566 contigs in final assembly and 500,159 of them were annotated to a total of 8,756 taxonomic labels


Across samples, we need to unify the taxonomic labels were used at colnames.


1. filtered contigs

There are a total of 105 plasmid contigs have hits with plasmid contigs (identity > 90, contig length > 800 bp, matched length/ contig length > 0.6) PLS_DB

By 2nd examination.

    35: chromosomem labeled as "chr" or "chr+pls"
    9: other plasmids, "plasmid" and no pTEF
    43: "pls+chr", 16 of which contains "pTEF"
      8 pTEF1: 1 pTEF1 + 7 pTEF1+VE14089
      6 pTEF2: 2 pTEF2; 2 pTEF2+ve; 2 ptef2+kb1
    17: "plasmid" + "pTEF"
      5: pTEF1
      8: ptef2
      3: pTEF3

    check for HIP:
    blast against HIP_pls contig: 4 passed filter
      k141_2180524: pls+chr
      k141_2242115: hits on Replicon, should be pTEF1
      k141_190777: pls+chr
      k141_237538: plasmid

To use, combine results from blast against pTEF/HIP

    pTEF1: 5 from plasmid and 2 from pls+chr
    pTEF2: 8 from plasmid and 1 from pls+chr
    pTEF3: 3 from plasmid
    HIP_pls: k141_237538, k141_190777, k141_2180524

    add 2 new contigs into "~/Documents/Work/project/HiC/CoAsm_150pe/Blast_result/plasmid_contig_checked.info"


18 contigs in total, linked with XXXX contigs, XXXX of which were annotated with XXX kraken labels

```{sh}

```


2. unfiltered contigs

  to use: 105 plasmid - "chr"
