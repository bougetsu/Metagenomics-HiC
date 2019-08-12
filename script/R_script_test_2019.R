setwd("C:/Work/HiC/1st/Result/links_count")
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
##########################
rc <- read.table("contig_stat", stringsAsFactors = F)

rc[,3] <-  as.numeric(gsub("read_count_", "", rc[,3]))

rc[,2] <-  as.numeric(gsub("length_", "", rc[,2]))
colnames(rc) <- c("name", "length", "count")
summary(rc[,3])
ggplot(rc, aes(x=count)) +
  geom_histogram(alha = 0.8, bins = 50) +
  xlim(c(100000, 4229541)) +
  geom_vline(aes(xintercept = 483458), color = "red")+
  geom_vline(aes(xintercept = 3229541), color = "blue")

####pad1
rc[rc[,1] == ">contig-100_569",]

####ptef2
rc[rc[,1] == ">contig-100_1705",]

###read data and transform###
#############################
###HiC links
links <- read.table("test_bam2links", header = F, stringsAsFactors = F)

str(links)
colnames(links) <- c("contig1", "contig2", "links")

links <- links %>%
    distinct(contig1, contig2, links)

inner_contig <- links[ which(links[,1] == links[,2]),]

###contig_stat
contig <- read.table("contig_stat", header = F, stringsAsFactors = F)

contig[,1] <- gsub(">", "", contig[,1])
contig[,2] <- as.numeric(gsub("length_", "", contig[,2]))
contig[,3] <- as.numeric(gsub("read_count_", "", contig[,3]))
colnames(contig) <- c("contig", "length", "count") 

###reads_on_contig

contig_reads <- read.table("test_bam2linksstat_reads_contigs", 
                           header = F, stringsAsFactors = F)
head(contig_reads)
colnames(contig_reads) <- c("contig", "Hic_reads")

contig_stat <- inner_join(contig, contig_reads, by = "contig") %>%
    filter(Hic_reads != 0)

####tax and cluster info
setwd("../20190218/")
tax <- read.table("bin_taxonomy.tab", sep = "\t", header = F, stringsAsFactors = F)
head(tax)
colnames(tax) <- c("cluster", "taxo")
tax[,1] <- gsub(".fa", "", tax[,1])
contig_bin <- read.table("metawrap_bins.contigs", sep = "\t",
                         header = F, stringsAsFactors = F)
head(contig_bin)
colnames(contig_bin) <- c("contig", "cluster")
cbintax <- right_join(tax, contig_bin, by = "cluster")

links_stat <- links %>%
    inner_join(contig_stat, by = c("contig1" = "contig")) %>%
    inner_join(contig_stat, by = c("contig2" = "contig")) %>%
    filter(length.x > 1000 & length.y > 1000) %>%
    mutate(prod_reads = count.x * count.y)


#################

l1 <- data.frame(c1= links[,1], c2 = links[,2])
l2 <- data.frame(c1= links[,2], c2 = links[,1])
idx.dup <- numeric()
for(i in 1:dim(l1)[1]){
    if(l1[i,1] != l1[i,2]){
        a <- which(l2[,1] == l1[i,1] & l2[,2] == l1[i,2])
        if(!is.na(a)){
            idx.dup <- c(idx.dup,a)
        }
    }
}


same_links <- links %>%
    filter(contig1 == contig2) %>%
    inner_join(contig_stat, by = c("contig1" = "contig")) %>%
    filter(length > 1000) %>%
    left_join(cbintax, by = c("contig1" = "contig"))
library(tidyr)
same_links <- same_links %>% replace_na(list(cluster = "unbinned", taxo = "unknown"))
same_links[,"taxo"] <- gsub("Bacteria;", "", same_links[,"taxo"])
################

same_links[,"taxo"] <- gsub("Actinobacteria;Actinobacteria;Bifidobacteriales;Bifidobacteriaceae", "Actinobacteria", same_links[,"taxo"])
same_links[,"taxo"] <- gsub("Bacteroidetes;Bacteroidia;Bacteroidales;Muribaculaceae", "Bacteroidetes", same_links[,"taxo"])
same_links[,"taxo"] <- gsub("Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecalis", "Firmicutes;E.faecalis", same_links[,"taxo"])
same_links[,"taxo"] <- gsub("Firmicutes;Clostridia;Clostridiales", "Firmicutes;Clostridia", same_links[,"taxo"])
same_links[,"taxo"] <- gsub("Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae;Faecalibaculum;Faecalibaculum rodentium", "Firmicutes;Erysipelotrichia", same_links[,"taxo"])
same_links[,"taxo"] <- gsub("Firmicutes;Erysipelotrichia;Erysipelotrichales;Erysipelotrichaceae", "Firmicutes;Erysipelotrichia", same_links[,"taxo"])
same_links[,"taxo"] <- gsub("Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae", "Bacteroidetes", same_links[,"taxo"])

ggplot(data = same_links, aes(x = count, y = links, color = taxo)) +
      geom_point(size = 1.5) +
    ylim(0, 2200000) +
    xlim(0,180000000) +
    theme_bw()

###############################
links_stat_tax <- links_stat %>%
    left_join(cbintax, by = c("contig1" = "contig")) %>%
    left_join(cbintax, by = c("contig2" = "contig")) %>%
    replace_na(list(cluster.x = "unbinned", taxo.x = "unknown",
               cluster.y = "unbinned", taxo.y = "unknown"))


links_ef <- links_stat_tax %>%
    filter(grepl("faecalis", taxo.x)) %>%
    mutate(type = ifelse(grepl("faecalis", taxo.y), "diff", "non_Ef")) %>%
    mutate(type=replace(type, contig1 == contig2, "same"))


ggplot(data = links_ef, aes(x = sqrt(prod_reads), y = links, color = type, shape = type)) +
    geom_point() +
    ylim(0, 10000)+
    theme_bw()

ggplot(data = links_ef[links_ef[,"type"] != "same",], aes(x = sqrt(prod_reads),
                                                          y = links, color = type, shape = type)) +
    geom_point() +
    theme_bw()




###try to find the abnormal high outer points
test <- arrange(links_ef[links_ef[,"type"] != "same",], desc(links)) %>%
    filter(type == "non_Ef")

test[1:30,]

################################################
#################################################


in2 <- links_stat %>%
#    filter(contig1 != contig2) %>%
    filter(length.x > 1000) %>%
    filter(length.y > 1000) %>%
    select(contig1, contig2, links, length.x, length.y)

############

stat <- cbintax %>%
    group_by(cluster) %>%
    summarise(n = n())

length(unique(in2[,1]))

cbintax %>%
    left_join(contig, by = "contig") %>%
    filter(length > 1000) %>%
    group_by(cluster) %>%
    summarise(n = n())


length(unique(c(in2[,1], in2[,2])))

#####length > 1000
sum(contig[,"length"] > 1000)
#7574

#####length > 1000 &HiC links
length(unique(c(in2[,1], in2[,2])))
#5368


####chech is contig > 1000bp with no HiC links
links_contig <- unique(links[,1])
length(links_contig)
### 26224

links_stat %>%
    filter(length.x > 1000) %>%
    distinct(contig1) %>%
    summarise(n = n())
##7503


links_stat %>%
    filter(length.x > 1000) %>%
    filter(length.y > 1000) %>%
    distinct(contig1) %>%
    summarise(n = n())
##7503


links_stat %>%
    filter(contig1 != contig2) %>%
    filter(length.x > 1000) %>%
    filter(length.y > 1000)%>%
    distinct(contig1) %>%
    summarise(n = n())
##5368

links_stat %>%
    filter(contig1 == contig2) %>%
    filter(length.x > 1000) %>%
    filter(length.y > 1000)%>%
    distinct(contig1) %>%
    summarise(n = n())
##7500



cbintax %>%
    left_join(contig, by = "contig") %>%
#    filter(length > 1000) %>%
    distinct(contig) %>%
    summarise(n = n())
### 2726


#####################


f3 <- in2 %>%
    left_join(cbintax, by = c("contig1" = "contig")) %>%
    left_join(cbintax, by = c("contig2" = "contig")) %>%
    replace_na(list(cluster.x = "unbinned", taxo.x = "unknown",
                    cluster.y = "unbinned", taxo.y = "unknown"))



stat2 <- f3 %>% select(contig1, cluster.x) %>%
    distinct(contig1, cluster.x) %>%
    group_by(cluster.x) %>%
    summarise(n= n())



# Check contigs bin/unbinned ----------------------------------------------


###################
f3 %>%
    filter(cluster.x != "unbinned") %>%  ###2762 binned
    filter(contig1 != contig2 ) %>%
    distinct(contig1) %>%
    summarise(n = n())
#2639    




f3 %>%
    filter(cluster.x != "unbinned") %>%  ###2762 binned
    filter(contig1 == contig2 ) %>%
    distinct(contig1) %>%
    summarise(n = n())
##2762


f3 %>%
    filter(cluster.x == "unbinned") %>%  ###4741 unbinned
    filter(contig1 != contig2 ) %>%  ##links with other contigs
    distinct(contig1) %>%
    summarise(n = n())
##2729


f3 %>%
    filter(cluster.x == "unbinned") %>%  ###4741 unbinned
    filter(contig1 == contig2 ) %>%  ##links with the self contigs
    distinct(contig1) %>%
    summarise(n = n())
##4738



#####binned contigs with other links
contig_bin_o <- f3 %>%
    filter(cluster.x != "unbinned") %>%  ###2762 binned
    filter(contig1 != contig2 ) %>%
    distinct(contig1)
contig_bin_o <- as.vector(contig_bin_o[,1])
length(contig_bin_o)
#2639
#####binned contigs with self links
contig_bin_s <- f3 %>%
    filter(cluster.x != "unbinned") %>%  ###2762 binned
    filter(contig1 == contig2 ) %>%
    distinct(contig1)
contig_bin_s <- as.vector(contig_bin_s[,1])
length(contig_bin_s)
#2672


ggplot(data = contig_stat, aes(x = count, y = Hic_reads)) +
    geom_point()

sum(contig[,"count"])


bob <- read.table("../201903/contig.blobplot", stringsAsFactors = F, sep = "\t", header = T)

dfb <- contig_stat %>%
    filter(contig %in% in2[,1] ) %>%
    left_join(bob, by = c("contig" = "seqid"))
library(skimr)
stats_bins <- dfb %>%
    select(length, count, cov_BLOBOLOGY.Mouse_Fecal_ShortGun.bowtie2.bam, Hic_reads, binned_yes_no) %>%
    group_by(binned_yes_no) %>% skim()

write.table(stats_bins, "../201903/stats_binsornot.tsv", sep = ",", quote = F, row.names = F)



dfb %>% group_by(binned_yes_no) %>%
    summarise(mean_len = mean(len), mean_cov = mean(cov_BLOBOLOGY.Mouse_Fecal_ShortGun.bowtie2.bam))


dfc <- dfb %>%
    select(length, count, cov_BLOBOLOGY.Mouse_Fecal_ShortGun.bowtie2.bam, Hic_reads, binned_yes_no) %>%
    gather(key = "feature", value = "value",
           length, count, cov_BLOBOLOGY.Mouse_Fecal_ShortGun.bowtie2.bam, Hic_reads)


ggplot(dfc, aes(x = log(value),  fill = binned_yes_no )) +
    geom_histogram() +
    facet_grid(binned_yes_no ~ feature, scales="free_x") +
    theme_bw()

ggplot(data = dfb, aes(x = gc, y = count/len, color = binned_yes_no)) +
    geom_point() +
    theme_bw()

ggplot(data = dfb, aes(x = gc, y = cov_BLOBOLOGY.Mouse_Fecal_ShortGun.bowtie2.bam, color = binned_yes_no)) +
    geom_point() +
    ylim(0, 10000) +
    theme_bw()


ggplot(data = dfb, aes(x = gc, y = len, color = binned_yes_no)) +
    geom_point()


ggplot(data = dfb, aes(x= log(len), fill = binned_yes_no)) +
    geom_histogram() +
    theme_bw()

ggplot(data = dfb, aes(x= log(count/len), fill = binned_yes_no)) +
    geom_histogram() +
    theme_bw()

ggplot(data = dfb, aes(x= log(cov_BLOBOLOGY.Mouse_Fecal_ShortGun.bowtie2.bam), fill = binned_yes_no)) +
    geom_histogram()


#########length > 1000 but interacted contigs < 1000
a <- setdiff(dfb[,1], in2[,1])




# HiC links with bins -----------------------------------------------------


#################


f4 <- f3 %>%
    filter(contig1 != contig2) %>%
    group_by(contig1, cluster.y, cluster.x)%>%
    summarise(HiC_links = sum(links), bin2_length = sum(length.y), n = n()) %>%
    mutate(ldl = HiC_links/bin2_length) %>%
    group_by(contig1) %>%
    arrange(desc(ldl), .by_group = TRUE)


f5 <- f4 %>%
    group_by(contig1) %>% top_n(2, ldl) %>%
    select(contig = contig1, cluster = cluster.y)


write.csv(f3, file = "links_contigs_tax.csv", quote = F, row.names = F)
write.csv(f4, file = "links_bins_tax.csv", quote = F, row.names = F)
#############
f3.2 <- in2 %>%
    left_join(f5, by = c("contig1" = "contig")) %>%
    left_join(f5, by = c("contig2" = "contig")) %>%
    replace_na(list(cluster.x = "unbinned", taxo.x = "unknown",
                    cluster.y = "unbinned", taxo.y = "unknown"))


f4.2 <- f3.2 %>%
    group_by(contig1, cluster.y, cluster.x)%>%
    summarise(HiC_links = sum(links), bin2_length = sum(length.y), n = n()) %>%
    mutate(ldl = HiC_links/bin2_length) %>%
    group_by(contig1) %>%
    arrange(desc(ldl), .by_group = TRUE)


f4.2 %>%
    group_by(contig1) %>% top_n(1, ldl)

f3.2 %>%
    filter(contig1 == "contig-100_1820")



f5.2 <- test12 <- test11 %>%
    group_by(contig1) %>% top_n(1, HiC_links)
test13 <- test12 %>% filter(cluster.y != cluster.x)

a <- paste(test11[,1],test11[,2])
b <- paste(test11[,2],test11[,1])
sum(a %in% b)



##############CHECK plasmid contig links###########


contig.pls <- c("569", "807", "2758")
contig.pls <- paste0("contig-100_", contig.pls)

links_pls <- links_stat_tax %>%
    filter(length.x >= 1000, length.y >= 1000) %>%
    filter(contig1 != contig2) %>%
    filter(contig1 %in% contig.pls)


l1 <- data.frame(c1= links_pls[,1], c2 = links_pls[,2], stringsAsFactors =F)
l2 <- data.frame(c1= links_pls[,2], c2 = links_pls[,1], stringsAsFactors =F)
idx.dup <- numeric()
for(i in 1:dim(l1)[1]){
    if(l1[i,1] != l1[i,2]){
        a <- which(l2[,1] == l1[i,1] & l2[,2] == l1[i,2])
        if(length(a) > 0){
            idx.dup <- c(idx.dup,a)
        }
    }
}

links_pls <- links_pls[-c(161, 270),]


tb_links_pls <- links_pls %>%
    group_by(contig1, cluster.y, cluster.x) %>%
    summarise(HiC_links = sum(links), bin2_length = sum(length.y), n = n(), ratio = HiC_links/bin2_length)
tb_links_pls
tax %>% left_join(tb_links_pls,by = c("cluster" = "cluster.y"))
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

library(pheatmap)
library("circlize")

adjacencyList2Matrix(pls_tb)
pheatmap(t(log2(adjacencyList2Matrix(pls_tb)+1)),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", 
         annotation_legend = T, legend = T)


tax$taxo <- gsub(".*;(*)", "\\1", tax$taxo)


write.table(tb_links_pls, file = "tb_links_plasmid", quote = F, row.names = F, sep = "\t")




###plot plasmid network
install.packages("GGally")
library(GGally)
install.packages("igraph")
library(igraph)


pls_cluster <- unique(pull(tb_links_pls,cluster.y))

edges <- links_stat_tax %>%
    filter(cluster.x %in% pls_cluster | contig1 %in% contig.pls) %>%
    filter(cluster.y %in% pls_cluster | contig2 %in% contig.pls) %>%
    filter(contig1 != contig2) %>%
    filter(cluster.x != "unbinned" | cluster.y != "unbinned") %>%
    filter

pls_contig <- unique(c(pull(edges,contig1), pull(edges,contig2)))


node <- contig_stat %>%
    filter(length > 1000) %>%
    left_join(cbintax, by = "contig") %>%
    replace_na(list(cluster = "unbinned", taxo = "unknown")) %>%
    mutate(plasmid = ifelse(contig %in% contig.pls, gsub("contig-100_","",contig), NA )) %>%
    filter(contig %in% pls_contig)


dim(node)

library(igraph)


edges[edges$contig1 == "contig-100_2671",]
edges[edges$contig2 == "contig-100_2671",]
edges <- edges[,1:3]
net <- graph_from_data_frame(d=edges, vertices=node, directed=F) 


colrs <- c( "#2166AC", "#92C5DE", "#F4A582", "#B2182B", "grey50")


t_col <- function(color, percent = 50, name = NULL) {
    #	  color = color name
    #	percent = % transparency
    #	   name = an optional name for the color
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100-percent)*255/100,
                 names = name)
    ## Save the color
    invisible(t.col)
    
}

colrs <- t_col(colrs)

V(net)$color <- colrs[factor(node$cluster)]
V(net)$frame.color <- ifelse(is.na(node$plasmid),NA, "black")
V(net)$size <- 1

net <- simplify(net, remove.multiple = T, remove.loops = T, edge.attr.comb=c(links="mean")) 
sum(degree(net) == 0)
l <- layout_with_lgl(net, root = "contig-100_569")
plot(net, vertex.label= node$plasmid,vertex.size = 4, vertex.label.cex = 0.6,edge.width = 0.5)

###############Same way check bin3C bin#################

###############
test11 <- contig_stat %>%
    filter(length > 1000) %>%
    left_join(cbintax, by = "contig") %>%
    replace_na(list(cluster = "unbinned", taxo = "unknown")) %>%
    mutate(plasmid = ifelse(contig %in% contig.pls, gsub("contig-100_","",contig), NA ))

test12 <- test11 %>%
    group_by(cluster, taxo) %>%
    summarise(abundance = sum(count))
write.table(test12, "piechart", sep = "\t", quote = F, row.names = F, col.names = F)



#########################











####################find bin3 and T11.OG1RF

setwd("../201903/")
T11 <- read.table("idba_1000@T11.tsv", sep = "\t", stringsAsFactors = F, header = F)
T11_contig <- unique(T11[,1])
OG1 <- read.table("idba_1000@OG1RF.tsv", sep = "\t", stringsAsFactors = F, header = F)
OG1_contig <- unique(OG1[,1])


length(T11_contig)
#281
length(OG1_contig)
#305

length(intersect(OG1_contig, T11_contig))
#273

length(setdiff(OG1_contig, T11_contig))
#32

length(setdiff(T11_contig, OG1_contig))
#8


T11[T11[,1] == "contig-100_569",]


##T11 plasmid
##NZ_GG688649.1

T11_p <- T11[T11[,2] == "NZ_GG688649.1",1]

dim(cbintax)
head(cbintax)


T11_d <- data.frame(contig = T11_contig, stringsAsFactors = F) %>%
    left_join(cbintax)

dTO <- data.frame(contig = c(T11_contig, OG1_contig), stringsAsFactors = F) %>%
    distinct(contig) %>%
    mutate(T11 = ifelse(contig %in% T11_contig,"Y",""), OG1RF = ifelse(contig %in% OG1_contig,"Y","")) %>%
    left_join(cbintax) %>%
    replace_na(list(cluster = "unbinned", taxo = "")) %>%
    left_join(contig_stat, by = "contig") %>%
    arrange(cluster, T11, OG1RF)

dTO[dTO[,"contig"] %in% T11_p,"T11"] <- "Plasmid"


dim(dTO)


write.csv(dTO, "contig_T11_og1.csv", quote = F, row.names = F)


dTO %>% filter(cluster == "bin.3") %>%
    filter(T11 == "" & OG1RF != "") %>%
    summarise(n = n())



head(cbintax)
dT1 <- cbintax %>% filter(cluster == "bin.3") %>%
    mutate(T11 = ifelse(contig %in% T11_contig,"Y",""), OG1RF = ifelse(contig %in% OG1_contig,"Y",""))%>%
    arrange(T11, OG1RF)
write.csv(dT1, "bin3_T11_og1.csv", quote = F, row.names = F)
########################################################

# pAD1 --------------------------------------------------------------------
pls_info <- read.csv("2019_03_05.tsv", stringsAsFactors = F, sep = "\t", header = T)
pls_info[,"Description_NUCCORE"] <- gsub(",", "", pls_info[,"Description_NUCCORE"])
cpls <- read.table("idba_800@PLS.tsv", stringsAsFactors = F)
colnames(cpls) <- c("contig", "ACC_NUCCORE","identity","align_length" , "qstart", "qend", "estart", "eend", "evalue", "bitscore",  "qq", "qv")
pls2 <- as.tibble(cpls) %>% left_join(pls_info, by = "ACC_NUCCORE") %>%
    select(contig, ACC_NUCCORE, identity,align_length, qstart, qend, estart, eend, evalue, Description_NUCCORE,taxon_name) %>%
    left_join(contig, by = "contig") %>%
    arrange(contig, qstart, identity)

write.csv(pls2, "contig_PLSDB.csv", quote = F, row.names = F)

pAD1 <- c("contig-100_569", "contig-100_807", "contig-100_2758")

pls3 <- pls2 %>%
    select(contig, taxon_name, length, count) %>%
    distinct(contig, length, count, Ef) %>%
    mutate(T11 = ifelse(contig %in% T11_contig,"Y",""),
           OG1RF = ifelse(contig %in% OG1_contig,"Y",""),
           pAD1 = ifelse(contig %in% pAD1,"Y","")) %>%
    left_join(cbintax)
write.csv(pls3, "contig_EF.csv", quote = F, row.names = F)


# metaphlan ---------------------------------------------------------------

setwd("C:\\Work/HiC/result/")

abund <- read.table("merged_abundance_table_species.txt", sep = "\t", stringsAsFactors = F, header = T)
colnames(abund) <- gsub("_Palmer_SG_S[0-9]_profile", "", colnames(abund))
colnames(abund) <- gsub("^X", "", colnames(abund))
colnames(abund) <- gsub("([0-9])([A-Z])", "\\1d_\\2", colnames(abund))
