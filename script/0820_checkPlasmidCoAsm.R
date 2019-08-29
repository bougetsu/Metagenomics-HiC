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
  group_by(ref)  %>%
  filter(len_match/len_contig > 0.6)%>%
  distinct(contig) %>%
  summarise(n = n())



coasm_v583 %>%
  filter(len_match/len_contig > 0.6) %>%
  filter(ref != "V583_chr") %>%
  distinct(contig) %>%
  write.table("V583_plasmid_contig.list", quote = F, row.names = F, col.names = F)


coasm_HIP <- read.table("CoAsm_HIP_genome.tsv", header = F, stringsAsFactors = F,
                            col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref))%>%
  filter(len_contig > 800, len_match > 800)

coasm_HIP %>%
  #group_by(ref) %>%
  filter(len_match/len_contig > 0.6)%>%
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
coasm_rep %>% filter(len_match/len_ref > 0.6) %>%
  left_join(pls_v583, by = "contig") %>%
  left_join(pls_HIP, by = "contig")


ptef_hip <- coasm_rep %>% filter(len_match > 300) %>%
  left_join(pls_v583, by = "contig") %>%
  left_join(pls_HIP, by = "contig") %>%
  distinct(contig) %>%
  pull(contig)

coasm_HIP %>% filter(contig %in% ptef_hip)
coasm_v583 %>% filter(contig %in% ptef_hip)


dim(coasm_rep)
contig_rep <- coasm_rep %>% filter(len_match > 300) %>%
  distinct(contig) %>%
  pull(contig)


##################
pls_info <- read.table("/ms/11/cong/data/PLS_db/1903/2019_03_05.tsv", stringsAsFactors = F, sep = "\t", quote = "", comment.char= "", header = T)
pls_name <- pls_info %>% select(ACC_NUCCORE, Description_NUCCORE, taxon_species_name)

coasm_PLS <- read.table("CoAsm_plsdb.tsv", header = F, stringsAsFactors = F,
                        col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref)) %>%
  filter(len_contig > 800, len_match/len_contig > 0.6)
coasm_PLS 
coasm_PLS <- coasm_PLS %>% group_by(contig, ref) %>%
  mutate(overlap.qstart = ifelse(qstart < qend, min(qstart), max(qstart)), 
         overlap.qend = ifelse(qstart < qend, max(qend), min(qend)),
         overlap.len_match = overlap.qend - overlap.qstart + 1) %>%
  filter(len_contig > 800, overlap.len_match/len_contig > 0.8)
coasm_PLS %>% group_by(contig) %>%
  arrange(desc(identity), desc(len_match)) %>%
  slice(1) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) 



coasm_PLS %>%
  filter(contig %in% contig_rep)%>%  left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) 


test <- coasm_PLS %>% slice(1)
test[2, ] <- coasm_PLS %>% slice(1)
test[3, ] <- coasm_PLS %>% slice(2)
test[2, "qstart"] <- 10
test[2, "qend"] <- 8180
test %>%group_by(contig, ref) %>%
  mutate(overlap.start = min(qstart), overlap.eng = max(qend))



coasm_v583 %>%
  filter(len_match/len_contig > 0.6) %>%
  filter(ref != "V583_chr") %>%
  group_by(contig) %>%
  mutate(ref = ifelse(n() > 1, "V583_pls", ref)) %>%
  distinct(contig, ref) 
  


pls_v583 <- coasm_v583 %>%
  filter(len_match/len_contig > 0.6) %>%
  filter(ref != "V583_chr") %>%
  group_by(contig) %>%
  mutate(ref = ifelse(n() > 1, "V583_pls", ref)) %>%
  distinct(contig, ref) 

pls_HIP <- coasm_HIP %>%
  filter(ref == "HIP_pls") %>%
  filter(len_match/len_contig > 0.6)%>%
  distinct(contig, ref)

coasm_rep <- coasm_rep %>% filter(len_match > 300) %>%
  distinct(contig) %>%
  mutate(Rep = c("rep9",
                 "rep1",
                 "repUS11",
                 "rep1",
                 "rep9",
                 "rep9"),
         ref = c("pTEF1",
                 "pTEF1",
                 "pTEF3",
                 "HIP_pls",
                 "pTEF2",
                 "gb|GG692639.1|"))
test <- pls_v583 %>% bind_rows(pls_HIP) %>%
  mutate(class = "plasmid") %>%
  full_join(coasm_rep, by = "contig")
