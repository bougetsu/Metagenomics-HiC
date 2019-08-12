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

setwd("/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/blast_result/")

V583_0d_100_v583 <- read.table("0V583_90_800.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref))

V583_0d_100_v583 %>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())


V583_13d_100_v583 <- read.table("13V583_90_800.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))%>%
  mutate(ref = sub_name(ref))

V583_13d_100_v583 %>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())

c_pls <- V583_13d_100_v583 %>%
  filter(ref != "ref|NC_004668.1|") %>%
  group_by(ref) %>%
  distinct(contig) %>%
  pull(contig) 

sort(c_pls[duplicated(c_pls)])


V583_13d_100_v583 %>%
  group_by(ref) %>%
  distinct(contig) %>%
  arrange(contig) 


pls_V583_13d <- V583_13d_100_v583 %>%
  group_by(ref) %>%
  distinct(contig)

####replicon
V583_0d_100_rep <- read.table("0V583_rep.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
####replicon
V583_13d_100_rep <- read.table("13V583_rep.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
V583_13d_100_rep %>%
  left_join(pls_V583_13d, by = "contig")
###pls_db
V583_0d_100_pls <- read.table("0V583_plsdb.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
V583_0d_100_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig) 
  
V583_13d_100_pls <- read.table("13V583_plsdb.tsv", header = F, stringsAsFactors = F,
                                col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
V583_13d_100_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig, Description_NUCCORE) %>%
  left_join(pls_V583_13d, by = "contig")

########################################3
######################################
###############150BP
######################################
##############################################
setwd("/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/150PE/blast_result/")

V583_0d_150_v583 <- read.table("filtered/0V583_150PE_genome_800_90.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref))

V583_0d_150_v583 %>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())


V583_13d_150_v583 <- read.table("filtered/13V583_150PE_genome_800_90.tsv", header = F, stringsAsFactors = F,
                                col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))%>%
  mutate(ref = sub_name(ref))

V583_13d_150_v583 %>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())

c_pls <- V583_13d_150_v583 %>%
  filter(ref != "ref|NC_004668.1|") %>%
  group_by(ref) %>%
  distinct(contig) %>%
  pull(contig) 

length(unique(c_pls[duplicated(c_pls)]))


V583_13d_100_v583 %>%
  group_by(ref) %>%
  distinct(contig) %>%
  arrange(contig) 


pls_V583_13d <- V583_13d_150_v583 %>%
  group_by(ref) %>%
  distinct(contig)

####replicon
V583_0d_150_rep <- read.table("0V583_150PE_rep.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
####replicon
V583_13d_150_rep <- read.table("13V583_150PE_rep.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
V583_13d_150_rep %>%
  left_join(pls_V583_13d, by = "contig")
###pls_db
V583_0d_150_pls <- read.table("filtered/0V583_150PE_plsdb_800_90.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
V583_0d_150_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig) 

V583_13d_150_pls <- read.table("13V583_150PE_plsdb.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
V583_13d_150_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig, Description_NUCCORE) %>%
  left_join(pls_V583_13d, by = "contig")
V583_13d_150_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig)






setwd("/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/blast_result/")

HIP_0d_100_HIP <- read.table("0HIP11704_genome.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref))

HIP_0d_100_HIP  %>%
  filter(len_match >= 800) %>%
  distinct(contig) %>%
  summarise(n = n())


HIP_13d_100_HIP <- read.table("13HIP11704_genome_90_800.tsv", header = F, stringsAsFactors = F,
                                col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))%>%
  mutate(ref = sub_name(ref))

HIP_13d_100_HIP%>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())

HIP_13d_100_HIP%>%
  filter(ref != "HIP_pls") %>%
  distinct(contig) %>%
  summarise(n = n())


pls_hip_13d <- HIP_13d_100_HIP %>%
  group_by(ref) %>%
  distinct(contig) 

####replicon
HIP_0d_100_rep <- read.table("0HIP11704_rep.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
HIP_0d_100_rep 
####replicon
HIP_13d_100_rep <- read.table("13HIP11704_rep.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
HIP_13d_100_rep %>%
  left_join(pls_hip_13d, by = "contig")
###pls_db
HIP_0d_100_pls <- read.table("0HIP11704_plsdb.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
HIP_0d_100_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig) 

HIP_13d_100_pls <- read.table("13HIP11704_plsdb.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
HIP_13d_100_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig)

########################################3
######################################
###############150BP
######################################
##############################################
setwd("/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/150PE/blast_result/")

HIP_0d_150_HIP <- read.table("filtered/0HIP11704_150PE_genome_800_90.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref))

HIP_0d_150_HIP%>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())


HIP_13d_150_HIP <- read.table("filtered/13HIP11704_150PE_genome_800_90.tsv", header = F, stringsAsFactors = F,
                                col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))%>%
  mutate(ref = sub_name(ref))

HIP_13d_150_HIP %>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())

HIP_13d_150_HIP%>%
  filter(ref != "HIP_pls") %>%
  distinct(contig) %>%
  summarise(n = n())


pls_HIP_13d <- HIP_13d_150_HIP %>%
  group_by(ref) %>%
  distinct(contig)

####replicon
HIP_0d_150_rep <- read.table("0HIP11704_150PE_rep.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
HIP_0d_150_rep
####replicon
HIP_13d_150_rep <- read.table("13HIP11704_150PE_rep.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
HIP_13d_150_rep %>%
  left_join(pls_HIP_13d, by = "contig")
###pls_db
HIP_0d_150_pls <- read.table("filtered/0HIP11704_150PE_plsdb_800_90.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
HIP_0d_150_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig) 

HIP_13d_150_pls <- read.table("13HIP11704_150PE_plsdb.tsv", header = F, stringsAsFactors = F,
                               col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
V583_13d_150_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig, Description_NUCCORE) %>%
  left_join(pls_V583_13d, by = "contig")
HIP_13d_150_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig)

# VE ----------------------------------------------------------------------

setwd("/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/blast_result/")

VE_0d_100_v583<- read.table("0VE14089_genome.tsv", header = F, stringsAsFactors = F,
                             col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref))

VE_0d_100_v583  %>%
  filter(len_match >= 800) %>%
  distinct(contig) %>%
  summarise(n = n())


VE_13d_100_v583 <- read.table("13VE14089_genome.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))%>%
  mutate(ref = sub_name(ref))

VE_13d_100_v583%>%
  filter(len_match >= 800) %>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())


pls_ve_13d <- VE_13d_100_v583 %>%
  filter(len_match >= 800)%>%
  group_by(ref) %>%
  distinct(contig) 

dup_contig <- pls_ve_13d[duplicated(pls_ve_13d[,1]),] %>% pull(contig)

pls_ve_13d %>% filter(contig %in% dup_contig)

VE_13d_100_v583 %>% filter(contig == "k119_43660")

####replicon
VE_0d_100_rep <- read.table("0VE14089_rep.tsv", header = F, stringsAsFactors = F,
                             col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
VE_0d_100_rep 
####replicon
VE_13d_100_rep <- read.table("13VE14089_rep.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
VE_13d_100_rep %>%
  left_join(pls_ve_13d, by = "contig")
###pls_db
VE_0d_100_pls <- read.table("0VE14089_plsdb.tsv", header = F, stringsAsFactors = F,
                             col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
VE_0d_100_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig) 

VE_13d_100_pls <- read.table("13VE14089_plsdb.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
VE_13d_100_pls %>% filter(len_match >= 800) %>% left_join(pls_name, by = c("ref"= "ACC_NUCCORE")) %>%
  distinct(contig)

########################################3
######################################
###############150BP
######################################
##############################################
setwd("/ms/11/cong/project/HiC/ASSEMBLY/MEGAHIT/150PE/blast_result/")

VE_0d_150_HIP <- read.table("filtered/0VE14089_150PE_HIP_800_90.tsv", header = F, stringsAsFactors = F,
                             col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref))

VE_0d_150_HIP%>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())


VE_0d_150_583<- read.table("filtered/0VE14089_150PE_V583_800_90.tsv", header = F, stringsAsFactors = F,
                            col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue")) %>%
  mutate(ref = sub_name(ref))

VE_0d_150_583%>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())

VE_13d_150_HIP <- read.table("filtered/13VE14089_150PE_HIP_800_90.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))%>%
  mutate(ref = sub_name(ref))

VE_13d_150_HIP %>%
#  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())


VE_13d_150_583 <- read.table("filtered/13VE14089_150PE_V583_800_90.tsv", header = F, stringsAsFactors = F,
                             col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))%>%
  mutate(ref = sub_name(ref))

VE_13d_150_583 %>%
  group_by(ref) %>%
  distinct(contig) %>%
  summarise(n = n())


pls_ve_13d <-VE_13d_150_583%>%
  group_by(ref) %>%
  distinct(contig)
VE_13d_150_583 %>% filter(contig == "k141_822148")
####replicon
VE_0d_150_rep <- read.table("0VE14089_150PE_rep.tsv", header = F, stringsAsFactors = F,
                             col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
VE_0d_150_rep
####replicon
VE_13d_150_rep <- read.table("13VE14089_150PE_rep.tsv", header = F, stringsAsFactors = F,
                              col.names = c("contig", "ref", "identity", "mis", "gap", "len_contig", "len_ref","len_match", "qstart", "qend", "sstart", "send", "evalue"))
VE_13d_150_rep %>%
  left_join(pls_ve_13d, by = "contig")
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
