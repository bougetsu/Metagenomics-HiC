if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")
library(Rsamtools)

#read in entire BAM file
setwd("/ms/11/cong/project/HiC/BINNING/150PE_coassembly/HiC_mapping/rep_res_mapping/")

bam <- scanBam("13HIP11704_hic2CoAsm2316_150pe.Rep.bam")

str(bam)

aln <- bam[[1]]
aln[1]
names(aln)
str(aln)

tab <- data.frame(reads = aln$qname, ref = aln$rname, pos = aln$pos, mate = aln$mrnm, mate.pos = aln$mpos, stringsAsFactors = F)

str(tab)
 tab %>% group_by(ref) %>%
   summarise(n = n()) %>%
   left_join(plasmid_contig_checked, by = c("ref" = "contig"))
 
 
rep9_ptef1 <- tab %>% filter(ref == "k141_813771") %>% pull(reads)

#rep9 on k141_813771: 3793    9099
tab %>% filter(ref == "k141_813771") 


which(aln$qname == rep9_ptef1)
aln$seq[1]
##got hit on VE14089   CHR
##Enterococcus faecalis strain VE14089 chromosome, complete genome	176	775	100%	2e-40	95.50%	CP039296.1
##NOT HIT WITH HIP plasmid, out of replicon range
## leave it alone

#rep1 on k141_2242115: 102     2174
rep1_ptef1 <- tab %>% filter(ref == "k141_2242115") %>% pull(reads)

reads_2242115 <- aln$seq[aln$qname %in% rep1_ptef1]
names(reads_2242115) <- aln$qname[aln$qname %in% rep1_ptef1]

which(aln$qname %in% rep1_ptef1)
aln$qname[aln$qname %in% rep1_ptef1]


writeXStringSet(reads_2242115, file = "13HIP_reads_2242115.fasta")

data.frame(reads = aln$qname[aln$qname %in% rep1_ptef1], seq = )


###K00364:132:H2VWMBBXY:4:1120:2696:40631 k141_2242115    3  k141_185471     1991
###Enterococcus faecalis strain VE14089 chromosome, complete genome	193	750	100%	1e-45	100.00%	CP039296.1
###K00364:132:H2VWMBBXY:4:1206:18802:18493 k141_2242115 2012 k141_2242115     2012
##################################
##LET CHECK VE STRAIN
####################################