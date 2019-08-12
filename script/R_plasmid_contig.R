library("dplyr")

setwd("/ms/11/cong/data/PLS_db/1903/")

info <- read.csv("2019_03_05.tsv", stringsAsFactors = F, sep = "\t", header = T)

head(info)

setwd("/ms/11/cong/project/HiC/data/raw_data/1st_mouse/ASSEMBLY/")

dat <- read.table("idbaud98.tsv", header = F, stringsAsFactors = F, sep = "\t")

head(dat)
colnames(dat) <- c("contig", "ACC_NUCCORE", "qstart", "qend", "estart", "eend", "evalue", "bitscore", "identity", "qq", "qv")

dat2 <- dat %>% left_join(info, by = "ACC_NUCCORE") %>%
  select(contig, ACC_NUCCORE, qstart, qend, estart, eend, evalue, identity, Description_NUCCORE, TaxonID_NUCCORE, TaxonID_NUCCORE, taxon_name)


library("dplyr")

setwd("/ms/11/cong/data/PLS_db/1903/")

info <- read.csv("2019_03_05.tsv", stringsAsFactors = F, sep = "\t", header = T)

head(info)

setwd("/ms/11/cong/project/HiC/data/raw_data/1st_mouse/")
cpls <- read.table("idba_800@PLS.tsv", stringsAsFactors = F)
head(cpls)

colnames(cpls) <- c("contig", "ACC_NUCCORE","identity","align_length" , "qstart", "qend", "estart", "eend", "evalue", "bitscore",  "qq", "qv")
pls2 <- cpls %>% left_join(info, by = "ACC_NUCCORE") %>%
  select(ACC_NUCCORE, identity,align_length, qstart, qend, estart, eend, evalue,  Description_NUCCORE, TaxonID_NUCCORE, TaxonID_NUCCORE, taxon_name) %>%
  arrange(taxon_name, contig, qstart, identity, align_length)