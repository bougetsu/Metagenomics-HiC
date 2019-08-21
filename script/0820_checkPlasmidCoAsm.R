###################
##check plasmid####
###################


library(dplyr)

V583_chr = "NC_004668.1"
pTEF1 = "NC_004669.1"
pTEF2 = "NC_004671.1"
pTEF3 = "NC_004669.1"
HIP_pls = "GG692659.1"

pls_info <- read.table("/ms/11/cong/data/PLS_db/1903/2019_03_05.tsv", stringsAsFactors = F, sep = "\t", quote = "", comment.char= "", header = T)
pls_name <- pls_info %>% select(ACC_NUCCORE, Description_NUCCORE)

sub_name <- Vectorize(function(x){
  if(x == "ref|NC_004668.1|") {
    return("V583_chr")
  }else if(x == "ref|NC_004669.1|"){
    return("pTEF1")
  }else if(x == "ref|NC_004671.1|"){
    return("pTEF2")
  }else if(x == "ref|NC_004670.1|"){
    return("pTEF3")
  }else if(x == "gb|GG692659.1|"){
    return("HIP_pls")
  }else{
    return(x)
  }
})

setwd("/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/CoAsm_150PE/blast_result/")

coasm_v583 <- read.table("CoAsm_V583_genome.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref)) %>%
  filter(len_contig > 800, len_match > 800)

coasm_v583 %>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())



coasm_HIP <- read.table("CoAsm_HIP_genome.tsv", header = F, stringsAsFactors = F,
                            col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref))%>%
  filter(len_contig > 800, len_match > 800)

coasm_HIP %>%
  #group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())

pls_v583 <- coasm_v583%>%
  group_by(ref) %>%
  distinct(contig)
pls_HIP <- coasm_HIP%>%
  group_by(ref) %>%
  distinct(contig)
coasm_HIP %>%
  group_by(ref) %>%
  filter(len_match/len_contig > 0.6)


####replicon
coasm_rep <- read.table("CoAsm_rep.tsv", header = F, stringsAsFactors = F,
                            col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
coasm_rep %>% filter(len_match > 300) %>%
  left_join(pls_v583, by = "contig") %>%
  left_join(pls_HIP, by = "contig")


ptef_hip <- coasm_rep %>% filter(len_match > 300) %>%
  left_join(pls_v583, by = "contig") %>%
  left_join(pls_HIP, by = "contig") %>%
  distinct(contig) %>%
  pull(contig)

coasm_HIP %>% filter(contig %in% ptef_hip)
coasm_v583 %>% filter(contig %in% ptef_hip)

####replicon
VE_13d_150_rep <- read.table("13VE14089_150PE_rep.tsv", header = F, stringsAsFactors = F,
                             col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
VE_13d_150_rep %>%
  left_join(pls_ve_13d, by = "contig") %>%
  left_join()
###pls_db
VE_0d_150_pls <- read.table("filtered/0VE14089_150PE_plsdb_800_90.tsv", header = F, stringsAsFactors = F,
                            col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
VE_0d_150_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig) 

VE_13d_150_pls <- read.table("13VE14089_150PE_plsdb.tsv", header = F, stringsAsFactors = F,
                             col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
V583_13d_150_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig, Description_NUCCORE) %>%
  left_join(pls_V583_13d, by = "contig")
VE_13d_150_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig)
