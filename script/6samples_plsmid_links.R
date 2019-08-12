#####
#find pls contif
####

library(dplyr)
setwd("/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/blast_result/")
day0 <- read.table("0V583_90_800.tsv", sep= "\t", stringsAsFactors = F)
colnames(day0) <-c("qseqid", "sseqid", "pident", "mismatch", "gapopen", "qlen", "slen", "length", "qstart", "qend", "sstart", "send", "evalue")
day0$sseqid <- gsub("ref\\|NC_004668.1\\|", "chr", day0$sseqid)
day0$sseqid <- gsub("ref\\|NC_004669.1\\|", "pTEF1", day0$sseqid)
day0$sseqid <- gsub("ref\\|NC_004670.1\\|", "pTEF3", day0$sseqid)
day0$sseqid <- gsub("ref\\|NC_004671.1\\|", "pTEF2", day0$sseqid)
day0  %>%
  group_by(sseqid) %>%
  summarize(n=n())



day13 <- read.table("13V583_90_800.tsv", sep= "\t", stringsAsFactors = F)
colnames(day13) <-c("qseqid", "sseqid", "pident", "mismatch", "gapopen", "qlen", "slen", "length", "qstart", "qend", "sstart", "send", "evalue")
day13$sseqid <- gsub("ref\\|NC_004668.1\\|", "chr", day13$sseqid)
day13$sseqid <- gsub("ref\\|NC_004669.1\\|", "pTEF1", day13$sseqid)
day13$sseqid <- gsub("ref\\|NC_004670.1\\|", "pTEF3", day13$sseqid)
day13$sseqid <- gsub("ref\\|NC_004671.1\\|", "pTEF2", day13$sseqid)
day13  %>%
  group_by(sseqid) %>%
  summarize(n=n())


library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)

setwd("/ms/11/cong/project/HiC/INITIAL_BINNING_0405/")
setwd("13V583/")

#############################
###read data and transform###
#############################
###HiC links
links <- read.table("13V583_hic_counts.tsv", header = F, stringsAsFactors = F)

str(links)
colnames(links) <- c("contig1", "contig2", "links")
links <- links %>%
  distinct(contig1, contig2, links)


####tax and cluster info
tax <- read.table("BIN_REFINEMENT/bin_taxonomy.tab", sep = "\t", header = F, stringsAsFactors = F)
head(tax)
colnames(tax) <- c("cluster", "taxo")
tax[,1] <- gsub(".fa", "", tax[,1])
contig_bin <- read.table("BIN_REFINEMENT/metawrap_50_10_bins.contigs", sep = "\t",
                         header = F, stringsAsFactors = F)
head(contig_bin)
colnames(contig_bin) <- c("contig", "cluster")
cbintax <- right_join(tax, contig_bin, by = "cluster")

links_tax <- links %>%
  left_join(cbintax, by = c("contig1" = "contig")) %>%
  left_join(cbintax, by = c("contig2" = "contig")) %>%
  filter(!is.na(cluster.y))


contig.pls <- day13 %>%
  filter(sseqid != "chr") %>%
  select(qseqid, sseqid)   %>%
  mutate(sseqid=replace(sseqid, qseqid %in% shared, "shared"))%>%
  distinct(qseqid, sseqid)
duplicated(contig.pls$qseqid)
shared <- c("k119_101703", "k119_143434", "k119_330295", "k119_334204", "k119_335017", "k119_403514", "k119_500141", "k119_50698", "k119_59995" )

#non-unique values when setting 'row.names': ‘k119_101703’, ‘k119_143434’, ‘k119_330295’, ‘k119_334204’, ‘k119_335017’, ‘k119_403514’, ‘k119_500141’, ‘k119_50698’, ‘k119_59995’ 


contig.pls_v <- contig.pls %>% pull(qseqid)

links_pls <- links_tax %>%
  filter(contig1 != contig2) %>%
  filter(contig1 %in% contig.pls_v)

tb_links_pls <- links_pls %>%
  group_by(contig1, cluster.y, cluster.x) %>%
  summarise(HiC_links = sum(links), n = n())
  #summarise(HiC_links = sum(links), bin2_length = sum(length.y), n = n(), ratio = HiC_links/bin2_length)
tb_links_pls
tax$taxo <- gsub(".*;(*)", "\\1", tax$taxo)
pls_tb<- tax %>% right_join(tb_links_pls,by = c("cluster" = "cluster.y")) %>%
  replace_na(list(cluster.x = "unbinned")) %>%
  select(contig1, taxo, HiC_links) %>%
  left_join(contig.pls, by = c("contig1" = "qseqid"))

library(pheatmap)
library("circlize")

mat <- adjacencyList2Matrix(pls_tb)
rownames(mat) == pls_tb$contig1
pls_row <- data.frame(plasmid = contig.pls$sseqid)
row.names(pls_row) <- contig.pls$qseqid

pheatmap(t(log2(adjacencyList2Matrix(pls_tb)+1)),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", 
         annotation_col = pls_row, legend = T)


########usining contig taxo


####tax and cluster info
contig_tax <- read.table("BIN_REFINEMENT/contig_taxonomy.tab", sep = "\t", header = F, stringsAsFactors = F)

colnames(contig_tax) <- c("contig", "taxo")
contig_tax <- contig_tax %>% filter(taxo != "") %>%
  mutate(taxo = gsub(".*;(*)", "\\1", taxo))
head(contig_tax)

links_contig_tax <- links %>%
  left_join(contig_tax, by = c("contig1" = "contig")) %>%
  left_join(contig_tax, by = c("contig2" = "contig")) %>%
  filter(!is.na(taxo.y))

head(links_contig_tax)

links_contig_pls <- links_contig_tax  %>%
  filter(contig1 != contig2) %>%
  filter(contig1 %in% contig.pls_v)
head(links_contig_pls)
tb_links_pls_contig <- links_contig_pls %>%
  group_by(contig1, taxo.y, taxo.x) %>%
  summarise(HiC_links = sum(links), n = n())
#summarise(HiC_links = sum(links), bin2_length = sum(length.y), n = n(), ratio = HiC_links/bin2_length)
tb_links_pls_contig

pls_tb_contig <- tb_links_pls_contig %>%
  left_join(contig.pls, by = c("contig1" = "qseqid")) %>%
  select(contig1, taxo.y, HiC_links)


mat2 <- adjacencyList2Matrix(pls_tb_contig)
rownames(mat2) == pls_tb_contig$contig1

pheatmap(t(log2(adjacencyList2Matrix(pls_tb_contig)+1)),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", 
         annotation_col = pls_row, legend = T,
         legend_breaks= c(0, 2, 4, 6, 8, 10, 12))
pheatmap(t(log2(adjacencyList2Matrix(pls_tb_contig)+1)),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         annotation_col = pls_row, legend = T)

t1 <- tb_links_pls %>%
  filter(contig1 == "contig-100_2758") %>%
  right_join(tax,by = c("cluster.y" = "cluster")) %>%
  replace_na(list(contig1 = "contig-100_2758", HiC_links = 0)) %>%
  select(contig1, taxo, HiC_links)
t2 <- tb_links_pls %>%
  filter(contig1 == "contig-100_569") %>%
  right_join(tax,by = c("cluster.y" = "cluster")) %>%
  replace_na(list(contig1 = "contig-100_569", HiC_links = 0)) %>%
  select(contig1, taxo, HiC_links)
t3 <- tb_links_pls %>%
  filter(contig1 == "contig-100_807") %>%
  right_join(tax,by = c("cluster.y" = "cluster")) %>%
  replace_na(list(contig1 = "contig-100_807", HiC_links = 0)) %>%
  select(contig1, taxo, HiC_links)
pls_tb <- t1 %>% bind_rows(t2) %>%
  bind_rows(t3)
install.packages("circlize")
library(pheatmap)
library("circlize")

adjacencyList2Matrix(pls_tb)
pheatmap(t(log(adjacencyList2Matrix(pls_tb)+1)),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", 
         annotation_legend = T, legend = T)


tax$taxo <- gsub(".*;(*)", "\\1", tax$taxo)
