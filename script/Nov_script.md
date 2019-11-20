1. manually modify E.f bin

    ```{sh}
    seqkit grep -f New_ef_v1_contig.list ../../final.contigs.fa >  New_ef_v1_contig.fa
    cat New_ef_v1_contig.fa manual.bin.19.fa > Manual.bin.19.v1.fa
    mkdir modified_ef_v1
    mv Manual.bin.19.v1.fa modified_ef_v1
    checkm lineage_wf -t 8 -x fa modified_ef_v1/ modified_ef_v1_checkm

    ```



| bin    | size | completedness | contamination rate | Strain heterogeneity |
|--------|------|---------------|--------------------|----------------------|
| before | 82   | 98.41         | 0.37               | 50                   |
| after  | 141  | 99.44         | 0.94               | 50.00                |
| after  | 338  | 99.44         | 1.69               | 66.67                |


only increase the contamination rate...

-------------------
*******************
-------------------


try those 6: from bin.141 and bin.199
k141_328065
k141_2292000
k141_914276
k141_372352
k141_1913480
k141_1659727
k141_1512295
k141_1015059

    ```{sh}
    mkdir modified_ef_v2
    seqkit grep -f New_ef_v2_contig.list ../../final.contigs.fa >  New_ef_v2_contig.fa
    cat New_ef_v2_contig.fa manual.bin.19.fa > modified_ef_v2/Manual.bin.19.v2.fa
    checkm lineage_wf -t 8 -x fa modified_ef_v2/ modified_ef_v2_checkm
    ```


| bin    | size | completedness | contamination rate | Strain heterogeneity |
|--------|------|---------------|--------------------|----------------------|
| before | 82   | 98.41         | 0.37               | 50                   |
| after  | 141  | 99.44         | 0.94               | 50.00                |
| after  | 149  | 99.44         | 0.94               | 50.00                |

We can keep those


................................
Next try

remove those from bin.141
k141_328065
k141_2292000
k141_914276
k141_372352
k141_1913480
k141_1659727
k141_1512295
```
/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/Taxify/metawrap/bin_141/remove_contig.v1


```
| bin    | size | completedness | contamination rate | Strain heterogeneity |
|--------|------|---------------|--------------------|----------------------|
| before | 744  | 75.66         | 3.21               | 33.33                |
| v1     | 737  | 75.66         | 3.21               | 33.33                |

k141_237718

| bin    | size | completedness | contamination rate | Strain heterogeneity |
|--------|------|---------------|--------------------|----------------------|
| before | 744  | 75.66         | 3.21               | 33.33                |
| v1     | 737  | 75.66         | 3.21               | 33.33                |
| v2     | 736  | 75.66         | 3.21               | 33.33                |

So its ok to move them to bin.19

..........................

Next

.........................

LETS REMOVE THE E.FACAELIS COMTAMINATION FIRST

bin.134

bin.134: Lactobacillus murinus

k141_1759014: Streptococcus

| bin    | size | completedness | contamination rate | Strain heterogeneity |
|--------|------|---------------|--------------------|----------------------|
| before | 744  | 99.48         | 0.00               | 0.00                 |
| v1     | 737  | 99.48         | 0.00               | 0.00                 |




2. Re-installed Anvio


```
conda activate anvio-master
anvi-activate-master


anvi-script-reformat-fasta ../../ASSEMBLY/MEGAHIT/CoAsm_150PE/final.contigs.fa -o contigs.fa --simplify-names
anvi-gen-contigs-database -f contigs.fa -o contigs.db -n 'CoAsm_150'

```
