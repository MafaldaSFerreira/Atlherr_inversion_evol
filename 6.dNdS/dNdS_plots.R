# Load packages ####
library(tidyverse)
library(biomaRt)

# Set workign directory ####
setwd("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Mutation_load/dNdS")

# Retrieve gene annotations from biomart ####
bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="charengus_gene_ensembl")
ensembl = useEnsembl(biomart="ensembl",dataset="charengus_gene_ensembl")
attr = c("ensembl_gene_id", "ensembl_transcript_id","external_gene_name", "chromosome_name","start_position",
         "end_position",
         "description")

# obtain annotations for all regions in the herring genome
regions <- getBM(attributes = attr, mart = bm)

colnames(regions)<-c("ENSEMBLID", "TRANSCRIPTID", "NAME", "Chrom", "Start", "End", "Description")

# Load input data ####
chr6_table <- read.table("data_2023-03-16/dnds.InvChr6.ML.txt", header=F)
colnames(chr6_table)<-c("transcript","dN_pop1","dS_pop1","dNdS_pop1","dN_pop2","dS_pop2","dNdS_pop2")
chr6_table$transcript<-str_remove(chr6_table$transcript, "/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/dNdS/CDS_seqs_by_transcript/alignments/")
chr6_table$ratio<-sqrt(chr6_table$dN_pop1)/sqrt(chr6_table$dN_pop2)
chr6_table$dNdSratio<-(chr6_table$dNdS_pop1)/(chr6_table$dNdS_pop2)

chr12_table <- read.table("data_2023-03-16/dnds.InvChr12.ML.txt", header=F)
colnames(chr12_table)<-c("transcript","dN_pop1","dS_pop1","dNdS_pop1","dN_pop2","dS_pop2","dNdS_pop2")
chr12_table$transcript<-str_remove(chr12_table$transcript, "/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/dNdS/CDS_seqs_by_transcript/alignments/")
chr12_table$ratio<-sqrt(chr12_table$dN_pop1)/sqrt(chr12_table$dN_pop2)
chr12_table$dNdSratio<-(chr12_table$dNdS_pop1)/(chr12_table$dNdS_pop2)

chr17_table <- read.table("data_2023-03-16/dnds.InvChr17.ML.txt", header=T)
colnames(chr17_table)<-c("transcript","dN_pop1","dS_pop1","dNdS_pop1","dN_pop2","dS_pop2","dNdS_pop2")
chr17_table$transcript<-str_remove(chr17_table$transcript, "/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/dNdS/CDS_seqs_by_transcript/alignments/")
chr17_table$ratio<-sqrt(chr17_table$dN_pop1)/sqrt(chr17_table$dN_pop2)
chr17_table$dNdSratio<-(chr17_table$dNdS_pop1)/(chr17_table$dNdS_pop2)

chr23_table <- read.table("data_2023-03-16/dnds.InvChr23.ML.txt", header=T)
colnames(chr23_table)<-c("transcript","dN_pop1","dS_pop1","dNdS_pop1","dN_pop2","dS_pop2","dNdS_pop2")
chr23_table$transcript<-str_remove(chr23_table$transcript, "/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/dNdS/CDS_seqs_by_transcript/alignments/")
chr23_table$ratio<-sqrt(chr23_table$dN_pop1)/sqrt(chr23_table$dN_pop2)
chr23_table$dNdSratio<-(chr23_table$dNdS_pop1)/(chr23_table$dNdS_pop2)

gw_table<-read.table("data_2023-03-16/dNdS_35i_61i_WG.txt", header=F)
colnames(gw_table)<-c("transcript","dN_pop1","dS_pop1","dNdS_pop1","dN_pop2","dS_pop2","dNdS_pop2")
gw_table$transcript<-str_remove(gw_table$transcript, "/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/dNdS/CDS_seqs_by_transcript/alignments/")
gw_table$ratio<-sqrt(gw_table$dN_pop1)/sqrt(gw_table$dN_pop2)
gw_table$dNdSratio<-(gw_table$dNdS_pop1)/(gw_table$dNdS_pop2)
gw_table<-merge(gw_table, regions, by.x="transcript", by.y="TRANSCRIPTID")

# Supplementary Figure 14 ####
# These are super extreme values that I think we should remove due to the possible problems with alignments...
gw_table<- gw_table %>% filter(dNdS_pop2 < 3)

chr6_table <- merge(chr6_table, regions, by.x="transcript", by.y="TRANSCRIPTID")
chr12_table <- merge(chr12_table, regions, by.x="transcript", by.y="TRANSCRIPTID")
chr17_table <- merge(chr17_table, regions, by.x="transcript", by.y="TRANSCRIPTID")
chr23_table <- merge(chr23_table, regions, by.x="transcript", by.y="TRANSCRIPTID")

chr6_table[2,]$NAME<-"ckap5"

chr6_dnds<-ggplot(chr6_table)+
  geom_point(aes(x=dNdS_pop1 , y=dNdS_pop2, color=abs(dNdSratio)!=1), alpha=0.5) +
  scale_color_manual(name="dNdSratio", values=setNames(c("orange","black"),c(T,F)))+
  geom_text(data=(chr6_table %>% filter(abs(dNdSratio)!=1)), aes(x=dNdS_pop1, y=dNdS_pop2+0.05, label=NAME), 
            hjust = "outward", size=2)+
  geom_abline(linetype="dashed", col="gray50")+
  xlab(expression(paste(d[N]/d[S]," Southern haplotype")))+
  ylab(expression(paste(d[N]/d[S]," Northern haplotype")))+
  ggtitle("Chromosome 6 inversion")+
  #geom_text(aes(x=0, y=0.90, label="38 genes"), hjust=0, size=5)+
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, size=7))

chr12_dnds<-ggplot(chr12_table)+
  geom_point(aes(x=dNdS_pop1 , y=dNdS_pop2, color=abs(dNdSratio)!=1), alpha=0.5) +
  scale_color_manual(name="dNdSratio", values=setNames(c("orange","black"),c(T,F)))+
  geom_text(data=(chr12_table %>% filter(abs(dNdSratio)!=1)), aes(x=dNdS_pop1, y=dNdS_pop2+0.05, label=NAME), 
            hjust = "outward", size=2)+
  geom_abline(linetype="dashed", col="gray50")+
  xlab(expression(paste(d[N]/d[S]," Southern haplotype")))+
  ylab(expression(paste(d[N]/d[S]," Northern haplotype")))+
  ggtitle("Chromosome 12 inversion")+
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=7),plot.title = element_text(hjust = 0.5, size=7))

chr17_dnds<-ggplot(chr17_table)+
  geom_point(aes(x=dNdS_pop1 , y=dNdS_pop2, color=abs(dNdSratio)!=1)) +
  scale_color_manual(name="dNdSratio", values=setNames(c("orange","black"),c(T,F)))+
  geom_abline(linetype="dashed", col="gray50")+
  xlab(expression(paste(d[N]/d[S]," Southern haplotype")))+
  ylab(expression(paste(d[N]/d[S]," Northern haplotype")))+
  ggtitle("Chromosome 17 inversion")+
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),plot.title = element_text(hjust = 0.5, size=7))

chr23_dnds<-ggplot(chr23_table)+
  geom_point(aes(x=dNdS_pop1 , y=dNdS_pop2, color=abs(dNdSratio)!=1)) +
  scale_color_manual(name="dNdSratio", values=setNames(c("orange","black"),c(T,F)))+
  geom_abline(linetype="dashed", col="gray50")+
  xlab(expression(paste(d[N]/d[S]," Southern haplotype")))+
  ylab(expression(paste(d[N]/d[S]," Northern haplotype")))+
  ggtitle("Chromosome 23 inversion")+
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=7),plot.title = element_text(hjust = 0.5, size=7)) 

all_plots<-grid.arrange(chr6_dnds, chr12_dnds, chr17_dnds, chr23_dnds, nrow=2)  

ggsave(all_plots, filename = "figures/dNdS_ratios_South_vs_North_inversions.pdf", height = 103, width=111, dpi=300, units="mm")


# Figure 6A. Boxplots ####

chr6_transcripts<-match(gw_table$transcript, chr6_table$transcript, nomatch=0) == 0
chr12_transcripts<-match(gw_table$transcript, chr12_table$transcript, nomatch=0) == 0
chr17_transcripts<-match(gw_table$transcript, chr17_table$transcript, nomatch=0) == 0
chr23_transcripts<-match(gw_table$transcript, chr23_table$transcript, nomatch=0) == 0

transcripts_to_remove <- chr6_transcripts + chr12_transcripts + chr17_transcripts +chr23_transcripts 
transcripts_to_keep <- transcripts_to_remove == 4 

gw_table_no_inversions<-gw_table[transcripts_to_keep,]

## Boxplots ####

chr6_table$region<-"chr6_inversion"
chr12_table$region<-"chr12_inversion"
chr17_table$region<-"chr17_inversion"
chr23_table$region<-"chr23_inversion"
gw_table_no_inversions$region<-"gw_table"

single_table_inversions<-rbind(chr6_table, chr12_table, chr17_table, chr23_table, gw_table_no_inversions)
single_table_inversions$Chrom<-as.character(single_table_inversions$Chrom)
single_table_inversions$region <- factor(single_table_inversions$region, levels=c("chr6_inversion", "chr12_inversion", "chr17_inversion", "chr23_inversion", "gw_table"))

# Remove the values pertaining to pop1 (which is calculated only with 35 individuals)
single_table_inversions[single_table_inversions$region=="gw_table",]$dNdS_pop1<-NA

dNdS_plot<-single_table_inversions %>%
  dplyr::select(transcript, region, dNdS_pop1, dNdS_pop2) %>%
  pivot_longer(cols=c(3,4), names_to="type" , values_to="dNdS") %>%
  ggplot()+
  geom_boxplot(aes(y=dNdS , x=region, fill=type), position="dodge", outlier.shape=NA)+
  scale_fill_manual(name="type", values=c("#D95F02","#1F78B4"))+
  ylim(0,1)+
  theme_classic() +
  theme(axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7))

ggsave(dNdS_plot, file="figures/dNdS_boxplot_North_vs_South_and_GW_20230808.pdf", height=50, width=198, units="mm")

## T-test ####
# I will do t-tests following Harringmeyer & Hoekstra. Within inversions, they used a two-sided t-test
# For genome-wide inversions, they calculated pNpS in sliding windows excluding the inversions and performed a one-sided t-test
# For now I only have box-plots for genome-wide using chr6 South and cht6 North, but we can collapse all 35 individuals to later obtain a background
chr6_ttest<-t.test(x = chr6_table$dNdS_pop1, y=chr6_table$dNdS_pop2, alternative = "two.sided")
chr12_ttest<-t.test(x = chr12_table$dNdS_pop1, y=chr12_table$dNdS_pop2, alternative = "two.sided")
chr17_ttest<-t.test(x = chr17_table$dNdS_pop1, y=chr17_table$dNdS_pop2, alternative = "two.sided")
chr23_ttest<-t.test(x = chr23_table$dNdS_pop1, y=chr23_table$dNdS_pop2, alternative = "two.sided")

ttest6_S_GW<-t.test(x = chr6_table$dNdS_pop1, y=gw_table_no_inversions$dNdS_pop1, alternative = "two.sided")
ttest6_N_GW<-t.test(x = chr6_table$dNdS_pop2, y=gw_table_no_inversions$dNdS_pop2, alternative = "two.sided")
ttest12_S_GW<-t.test(x = chr12_table$dNdS_pop1, y=gw_table_no_inversions$dNdS_pop1, alternative = "two.sided")
ttest12_N_GW<-t.test(x = chr12_table$dNdS_pop2, y=gw_table_no_inversions$dNdS_pop2, alternative = "two.sided")
ttest17_S_GW<-t.test(x = chr17_table$dNdS_pop1, y=gw_table_no_inversions$dNdS_pop1, alternative = "two.sided")
ttest17_N_GW<-t.test(x = chr17_table$dNdS_pop2, y=gw_table_no_inversions$dNdS_pop2, alternative = "two.sided")
ttest23_S_GW<-t.test(x = chr23_table$dNdS_pop1, y=gw_table_no_inversions$dNdS_pop1, alternative = "two.sided")
ttest23_N_GW<-t.test(x = chr23_table$dNdS_pop2, y=gw_table_no_inversions$dNdS_pop2, alternative = "two.sided")

t.test_dnds_results<-data.frame(comparisons=c(chr6_ttest$data.name,
                                              chr12_ttest$data.name,
                                              chr17_ttest$data.name,
                                              chr23_ttest$data.name,
                                              ttest6_S_GW$data.name,
                                              ttest6_N_GW$data.name,
                                              ttest12_S_GW$data.name,
                                              ttest12_N_GW$data.name,
                                              ttest17_S_GW$data.name,
                                              ttest17_N_GW$data.name,
                                              ttest23_S_GW$data.name,
                                              ttest23_N_GW$data.name),
                                tvalue=c(chr6_ttest$statistic,
                                         chr12_ttest$statistic,
                                         chr17_ttest$statistic,
                                         chr23_ttest$statistic,
                                         ttest6_S_GW$statistic,
                                         ttest6_N_GW$statistic,
                                         ttest12_S_GW$statistic,
                                         ttest12_N_GW$statistic,
                                         ttest17_S_GW$statistic,
                                         ttest17_N_GW$statistic,
                                         ttest23_S_GW$statistic,
                                         ttest23_N_GW$statistic),
                                df=c(chr6_ttest$parameter,
                                     chr12_ttest$parameter,
                                     chr17_ttest$parameter,
                                     chr23_ttest$parameter,
                                     ttest6_S_GW$parameter,
                                     ttest6_N_GW$parameter,
                                     ttest12_S_GW$parameter,
                                     ttest12_N_GW$parameter,
                                     ttest17_S_GW$parameter,
                                     ttest17_N_GW$parameter,
                                     ttest23_S_GW$parameter,
                                     ttest23_N_GW$parameter),
                                pvalue=c(chr6_ttest$p.value,
                                         chr12_ttest$p.value,
                                         chr17_ttest$p.value,
                                         chr23_ttest$p.value,
                                         ttest6_S_GW$p.value,
                                         ttest6_N_GW$p.value,
                                         ttest12_S_GW$p.value,
                                         ttest12_N_GW$p.value,
                                         ttest17_S_GW$p.value,
                                         ttest17_N_GW$p.value,
                                         ttest23_S_GW$p.value,
                                         ttest23_N_GW$p.value))

write.table(t.test_dnds_results, file="t_test_results/t_test_results_dnds_North_vs_South_vs_GW.txt", quote=F, sep="\t", col.names = T, row.names = F)
write.table(t.test_dnds_results, file="t_test_results/t_test_results_dnds_North_vs_South_vs_GW_20230808.txt", quote=F, sep="\t", col.names = T, row.names = F)

## write_outputs ####

bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="charengus_gene_ensembl")
ensembl = useEnsembl(biomart="ensembl",dataset="charengus_gene_ensembl")
attr = c("ensembl_gene_id", "ensembl_transcript_id","external_gene_name", "chromosome_name","start_position",
         "end_position",
         "description")

# obtain annotations for all regions in the herring genome
regions <- getBM(attributes = attr, mart = bm)

colnames(regions)<-c("ENSEMBLID", "TRANSCRIPTID", "NAME", "Chrom", "Start", "End", "Description")

chr6_table$chr <- "chr6"
chr12_table$chr <- "chr12"
chr17_table$chr <- "chr17"
chr23_table$chr <- "chr23"

output_table<-rbind(chr6_table, chr12_table, chr17_table, chr23_table)

output_table_annot<-merge(output_table, regions, by.x="transcript", by.y="TRANSCRIPTID")

write.table(output_table_annot, file ="dnds_results_all_inversions_2023-03-16.txt", col.names = T, row.names = F, quote = F, sep="\t")
