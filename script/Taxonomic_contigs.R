###############################
###taxonomic annotation of kraken, diamond_refpro, diamond_nt, metawrap

library(dplyr)
library(stringr)
setwd("/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/")


###################refpro

contig_id <- read.table("Taxify/refpro/blob_taxify.assembly.vs.refpro.1e25.diamond.taxified.out", stringsAsFactors = F)
colnames(contig_id) <- c("contig", "taxid", "EScore", "ssid")

taxid <- read.table("Taxify/refpro/taxid_name.txt", sep = "\t", stringsAsFactors = F, skip = 1, header = T)

ref_pro <- contig_id %>% left_join(taxid, by = "taxid") %>%
  select(contig, taxname, lineage)

###################kraken

kraken <- read.table("Taxify/Kraken/final.contigs.kraLabel", sep = "\t", stringsAsFactors = F, header = F)
colnames(kraken) <- c("contig", "kra_label")



#################bin_contig
bin <- read.table("Taxify/metawrap/bin_contig.tsv", sep = "\t", stringsAsFactors = F, header = T)

nonef <- bin %>% full_join(kraken, by = "contig") %>%
  full_join(ref_pro, by = "contig") %>%
  filter(cluster != "bin.19")

dim(nonef)
nonef <- nonef %>% filter(!kra_label == "cellular organisms;Bacteria") %>%
  filter(!kra_label == "cellular organisms")
dim(nonef)

nonef <- nonef %>% select(contig, cluster, kra_label, taxname)
idx_kef <- which(grepl("Enterococcus",nonef$kra_label))
nonef[1264,]
list1 <- nonef[idx_kef[1:30],]
nonef <- nonef[-(1:9000),]
dim(nonef)
idx_kef <- which(grepl("Enterococcus",nonef$kra_label))
list2 <- nonef[idx_kef,]
write.table(nonef, "Taxify/metawrap/kra_nonef.tsv", sep = "\t", row.names = F, quote = T, col.names = F)
write.table(list1, "Taxify/metawrap/kra_check_1", sep = "\t", row.names = F, quote = F)
check_kra <- rbind(list1, list2)



taxid <- read.table("Taxify/metawrap/taxid2name_nt_check", sep = "\t", stringsAsFactors = F,header = T, skip = 1)
colnames(taxid) <- c("taxid", "taxname")
nt_check_contigs <- read.table("Taxify/metawrap/nt_check_contigs", sep = "\t", stringsAsFactors = F)
nt_check_contigs <- nt_check_contigs[,c(1, 2, 3)]
colnames(nt_check_contigs) <- c("contig", "taxid", "score")

nt_check_contigs<- nt_check_contigs %>% left_join(taxid, by = "taxid") %>%
  select(contig, nt_label = taxname, taxid, score)

nt_check_contigs_top3 <- nt_check_contigs %>% group_by(contig) %>%
  arrange(desc(score)) %>% slice(1:3)

nt_check_contigs_top3  %>% filter(is.na(nt_label)) %>%
  pull(taxid)

nt_check_contigs %>% group_by(contig) %>%
  arrange(desc(score)) %>% slice(1) %>% right_join(check_kra, by = "contig") %>%
  select(contig, cluster, kra_label, nt_label, ref_pro) %>%
  write.table("Taxify/metawrap/To_check_ef.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

check_kra <- read.table("Taxify/metawrap/To_check_ef.tsv", sep = "\t", header = T)