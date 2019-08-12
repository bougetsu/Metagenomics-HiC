#################
##compare bin evaluation: metawrap bin_refinement, DAS, metabat2, maxbin2, bin3C
################


install.packages("doMC")
install.packages("data.table")
library(doMC)
library("data.table")
library(dplyr)
library(ggplot2)
install.packages("gridExtra")
library(gridExtra)

setwd("/ms/11/cong/project/HiC/INITIAL_BINNING_0405/")
for( s in dir()) {
  com_checkm(s)
}
setwd("/ms/11/cong/project/HiC/INITIAL_BINNING_CO_0426/")
com_checkm <- function(ss){
  setwd(ss)
  Metabat2 <- read.delim("refinenment_5k//metabat2_bins.stats", sep = "\t", header = TRUE) %>%
    select(bin, completeness, contamination, N50, GC, size) %>%
    rename(Completeness = completeness, Contamination = contamination)
  Maxbin2 <- read.delim("refinenment_5k/maxbin2_bins.stats", sep = "\t", header = TRUE) %>%
    select(bin, completeness, contamination, N50, GC, size) %>%
    rename(Completeness = completeness, Contamination = contamination)
  bin3C <- read.delim("refinenment_5k/5k_fasta.stats", sep = "\t", header = TRUE) %>%
    select(bin, completeness, contamination, N50, GC, size) %>%
    rename(Completeness = completeness, Contamination = contamination)
  Metawrap <- read.delim("refinenment_5k/metawrap_50_10_bins.stats", sep = "\t", header = TRUE) %>%
    select(bin, completeness, contamination, N50, GC, size) %>%
    rename(Completeness = completeness, Contamination = contamination)
  DASTool <- read.delim("DAS/checkM/das_checkm.results.tab", sep = "\t", header = TRUE) %>%
    select(Bin.Id, Completeness,Contamination, N50..contigs. , GC, Genome.size..bp.) %>%
    rename(bin = Bin.Id, size = Genome.size..bp., N50 = N50..contigs.)
  
  pbmethods <- list(Metabat2 = Metabat2, Maxbin2 = Maxbin2, bin3C = bin3C, Metawrap = Metawrap, DASTool = DASTool)
  
  methods <- names(pbmethods)
  for(i in 1:length(pbmethods)){
    pbmethods[[i]]$Method <- methods[i]
  }
  result_table <- do.call(rbind.data.frame, pbmethods)
  
  result_table$Method <- factor(result_table$Method, levels = c("Metabat2", "Maxbin2", "bin3C", "Metawrap", "DASTool"))
  tmp_wide <- result_table %>% group_by(Method) %>%
    filter(Contamination <= 10) %>%
    summarize(`>90%` = sum(Completeness>90), `>80%` = sum(Completeness>80 & Completeness<=90),
              `>70%` = sum(Completeness>70 & Completeness<=80), `>60%` = sum(Completeness>60 & Completeness<=70),
              `>50%` = sum(Completeness>50 & Completeness<=60), Sum = sum(Completeness>=50))
  
  plot_table <- melt(tmp_wide,id.vars = 'Method',
                     measure.vars = c('>90%','>80%','>70%', '>60%', '>50%'),
                     value.name = 'bin',variable.name ="Completeness")
  plot_table$title <- "CheckM Bin Statistics (<= 10% Contamination)"
  plot_table$Completeness <- factor(plot_table$Completeness,levels = rev(c('>90%','>80%','>70%', '>60%', '>50%')))
  colors <- rev(c("#08306B","#1664AB","#4A97C9","#93C4DE", "#eaf2fd"))
  
  
  tmp_wide2 <- result_table %>% group_by(Method) %>%
    filter(Contamination <= 5) %>%
    summarize(`>90%` = sum(Completeness>90), `>80%` = sum(Completeness>80 & Completeness<=90),
              `>70%` = sum(Completeness>70 & Completeness<=80), `>60%` = sum(Completeness>60 & Completeness<=70),
              `>50%` = sum(Completeness>50 & Completeness<=60), Sum = sum(Completeness>=50))
  
  plot_table2 <- melt(tmp_wide2,id.vars = 'Method',
                     measure.vars = c('>90%','>80%','>70%', '>60%', '>50%'),
                     value.name = 'bin',variable.name ="Completeness")
  plot_table2$title <- "CheckM Bin Statistics (<= 5% Contamination)"
  plot_table2$Completeness <- factor(plot_table2$Completeness,levels = rev(c('>90%','>80%','>70%', '>60%', '>50%')))
  
  
  tmp2 <- result_table %>% 
    filter(Contamination <= 10, Completeness>=50) %>%
    group_by(Method)  %>%
    arrange(desc(Completeness)) %>%
    mutate(Completeness_Rank = order(Completeness, decreasing = T))
  
  tmp3 <- result_table %>% 
    filter(Contamination <= 10, Completeness>=50) %>%
    group_by(Method)  %>%
    arrange(Method, Contamination) %>%
    mutate(Contamination_Rank = order(Method, Contamination, decreasing = F))
  
  tmp4 <- result_table %>% 
    filter(Contamination <= 10, Completeness>=50)
  par(mfrow = c(1, 3))
  p1 <- ggplot(plot_table, aes(Method, bin, fill=Completeness)) + geom_bar(stat="identity", position="stack")  +
    facet_grid( ~ title)  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = colors) +  scale_x_discrete(limits = methods)+
    ylim(0, 80)
  p2 <- ggplot(plot_table2, aes(Method, bin, fill=Completeness)) + geom_bar(stat="identity", position="stack")  +
    facet_grid( ~ title)  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = colors) +  scale_x_discrete(limits = methods)+
    ylim(0, 80)
  p3 <- ggplot(tmp2, aes(x = Completeness_Rank, y = Completeness, color = Method,group = Method))+
    geom_point() +
    geom_line(size = 1.2) +
    theme_bw()
  p4 <- ggplot(tmp3, aes(x = Contamination_Rank, y = Contamination, color = Method,group = Method))+
    geom_point() +
    geom_line(size = 1.2) +
    theme_bw()
  p5 <- ggplot(tmp4, aes(x = size, color = Method,fill=Method)) +
    geom_histogram(bins = 30,  alpha=0.2, position="identity")+
    facet_grid(Method ~ .)
  p6 <- ggplot(tmp4, aes(x = N50, color = Method,fill=Method)) +
    geom_histogram(bins = 30,  alpha=0.2, position="identity")+
    facet_grid(Method ~ .)
  setwd("..")
  pdf(paste0(ss, "_check_dastool_stats.pdf"), useDingbats=FALSE,width = 18, height = 12)
  grid.arrange(
    grobs = list(p1, p2, p3, p4, p5, p6),
    widths = c(1, 1, 1, 1),
    top = ss,
    layout_matrix = rbind(c(1, 2, 5, 6),
                          c(3, 4, 5, 6))
  )
  dev.off()
  write.csv(tmp_wide, file = paste0(ss, "_goodbin.csv"), quote = F, row.names = F)
}






com_checkm <- function(ss){
  setwd(ss)
  Metabat2 <- read.delim("refinenment_5k//metabat2_bins.stats", sep = "\t", header = TRUE) %>%
    select(bin, completeness, contamination, N50, GC, size) %>%
    rename(Completeness = completeness, Contamination = contamination)
  Maxbin2 <- read.delim("refinenment_5k/maxbin2_bins.stats", sep = "\t", header = TRUE) %>%
    select(bin, completeness, contamination, N50, GC, size) %>%
    rename(Completeness = completeness, Contamination = contamination)
  bin3C <- read.delim("refinenment_5k/5k_fasta.stats", sep = "\t", header = TRUE) %>%
    select(bin, completeness, contamination, N50, GC, size) %>%
    rename(Completeness = completeness, Contamination = contamination)
  Metawrap <- read.delim("refinenment_5k/metawrap_50_10_bins.stats", sep = "\t", header = TRUE) %>%
    select(bin, completeness, contamination, N50, GC, size) %>%
    rename(Completeness = completeness, Contamination = contamination)
  DASTool <- read.delim("DAS/checkM/das_checkm.results.tab", sep = "\t", header = TRUE) %>%
    select(Bin.Id, Completeness,Contamination, N50..contigs. , GC, Genome.size..bp.) %>%
    rename(bin = Bin.Id, size = Genome.size..bp., N50 = N50..contigs.)
  
  pbmethods <- list(Metabat2 = Metabat2, Maxbin2 = Maxbin2, bin3C = bin3C, Metawrap = Metawrap, DASTool = DASTool)
  
  methods <- names(pbmethods)
  for(i in 1:length(pbmethods)){
    pbmethods[[i]]$Method <- methods[i]
  }
  result_table <- do.call(rbind.data.frame, pbmethods)
  
  result_table$Method <- factor(result_table$Method, levels = c("Metabat2", "Maxbin2", "bin3C", "Metawrap", "DASTool"))
  tmp_wide <- result_table %>% group_by(Method) %>%
    filter(Contamination <= 10) %>%
    summarize(`>90%` = sum(Completeness>90), `>80%` = sum(Completeness>80 & Completeness<=90),
              `>70%` = sum(Completeness>70 & Completeness<=80), `>60%` = sum(Completeness>60 & Completeness<=70),
              `>50%` = sum(Completeness>50 & Completeness<=60), Sum = sum(Completeness>=50))
  
  plot_table <- melt(tmp_wide,id.vars = 'Method',
                     measure.vars = c('>90%','>80%','>70%', '>60%', '>50%'),
                     value.name = 'bin',variable.name ="Completeness")
  plot_table$title <- "CheckM Bin Statistics (<= 10% Contamination)"
  plot_table$Completeness <- factor(plot_table$Completeness,levels = rev(c('>90%','>80%','>70%', '>60%', '>50%')))
  colors <- rev(c("#08306B","#1664AB","#4A97C9","#93C4DE", "#eaf2fd"))
  
  
  tmp_wide2 <- result_table %>% group_by(Method) %>%
    filter(Contamination <= 5) %>%
    summarize(`>90%` = sum(Completeness>90), `>80%` = sum(Completeness>80 & Completeness<=90),
              `>70%` = sum(Completeness>70 & Completeness<=80), `>60%` = sum(Completeness>60 & Completeness<=70),
              `>50%` = sum(Completeness>50 & Completeness<=60), Sum = sum(Completeness>=50))
  
  plot_table2 <- melt(tmp_wide2,id.vars = 'Method',
                      measure.vars = c('>90%','>80%','>70%', '>60%', '>50%'),
                      value.name = 'bin',variable.name ="Completeness")
  plot_table2$title <- "CheckM Bin Statistics (<= 5% Contamination)"
  plot_table2$Completeness <- factor(plot_table2$Completeness,levels = rev(c('>90%','>80%','>70%', '>60%', '>50%')))
  
  
  tmp2 <- result_table %>% 
    filter(Contamination <= 10, Completeness>=50) %>%
    group_by(Method)  %>%
    arrange(desc(Completeness)) %>%
    mutate(Completeness_Rank = order(Completeness, decreasing = T))
  
  tmp3 <- result_table %>% 
    filter(Contamination <= 10, Completeness>=50) %>%
    group_by(Method)  %>%
    arrange(Method, Contamination) %>%
    mutate(Contamination_Rank = order(Method, Contamination, decreasing = F))
  
  tmp4 <- result_table %>% 
    filter(Contamination <= 10, Completeness>=50)
  par(mfrow = c(1, 3))
  p1 <- ggplot(plot_table, aes(Method, bin, fill=Completeness)) + geom_bar(stat="identity", position="stack")  +
    facet_grid( ~ title)  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = colors) +  scale_x_discrete(limits = methods)+
    ylim(0, 200)
  p2 <- ggplot(plot_table2, aes(Method, bin, fill=Completeness)) + geom_bar(stat="identity", position="stack")  +
    facet_grid( ~ title)  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = colors) +  scale_x_discrete(limits = methods)+
    ylim(0, 200)
  p3 <- ggplot(tmp2, aes(x = Completeness_Rank, y = Completeness, color = Method,group = Method))+
    geom_point() +
    geom_line(size = 1.2) +
    theme_bw()
  p4 <- ggplot(tmp3, aes(x = Contamination_Rank, y = Contamination, color = Method,group = Method))+
    geom_point() +
    geom_line(size = 1.2) +
    theme_bw()
  p5 <- ggplot(tmp4, aes(x = size, color = Method,fill=Method)) +
    geom_histogram(bins = 30,  alpha=0.2, position="identity")+
    facet_grid(Method ~ .)
  p6 <- ggplot(tmp4, aes(x = N50, color = Method,fill=Method)) +
    geom_histogram(bins = 30,  alpha=0.2, position="identity")+
    facet_grid(Method ~ .)
  setwd("..")
  pdf(paste0(ss, "_check_dastool_stats.pdf"), useDingbats=FALSE,width = 18, height = 12)
  grid.arrange(
    grobs = list(p1, p2, p3, p4, p5, p6),
    widths = c(1, 1, 1, 1),
    top = ss,
    layout_matrix = rbind(c(1, 2, 5, 6),
                          c(3, 4, 5, 6))
  )
  dev.off()
  write.csv(tmp_wide, file = paste0(ss, "_goodbin.csv"), quote = F, row.names = F)
}
