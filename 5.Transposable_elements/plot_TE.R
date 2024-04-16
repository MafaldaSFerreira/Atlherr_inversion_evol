# Plot TE masking results for inversions #
# Working directory ####
setwd("~/Documents/Postdoc/Project_Herring/Inversion_Project/RepeatMasker/find_inversions")

# Libraries ####
library(ggplot2)
#install.packages("rstatix")
library(rstatix)
library(tidyverse)
#install.packages("ggsignif")
library(ggsignif)

# Input file ####
table<-read.table("results_for_R.txt", header=T)

table$Database<-factor(table$Database,levels = c("HerringCombinedDB","Actinopterygii"))
table$Chromosome<-factor(table$Chromosome,levels = c("Genome","Chr6","Chr12","Chr17","Chr23"))
table$Genotype<-factor(table$Genotype,levels = c("S", "N", "NS"))

# Make the whole-genome a single data-point
table[table$Chromosome=="Genome",]$Genotype<-"N"

masked_fish<-table %>% 
  filter(Genotype!="NS") %>%
  filter(Database=="Actinopterygii") %>%
  ggplot()+
  geom_boxplot(aes(y=Masked_Proportion, x=Chromosome, color=Genotype))+
  ylim(0, 45)+
  theme_bw()+
  ylab("Masked Proportion (%)")+
  theme(axis.text.y = element_text(colour="grey10",size=15),
        axis.title.y = element_text(colour="black",size=15),
        axis.text.x = element_text(colour="grey10",size=15),
        axis.title.x = element_text(colour="black",size=15),
        title = element_text(size=15),
        legend.text=element_text(size=15))+
  ggtitle("Database: Actinopterygii")

masked_herring<-table %>% 
  filter(Genotype!="NS") %>%
  filter(Database=="HerringCombinedDB") %>%
  ggplot(aes(y=Masked_Proportion, x=Chromosome, fill=Genotype))+
  geom_boxplot()+
  scale_fill_manual(values=c("#2A2B78", "#FFC000")) +
  ylim(0, 45)+
  theme_bw()+
  ylab("Masked Proportion (%)")+
  theme(axis.text.y = element_text(colour="grey10",size=15),
        axis.title.y = element_text(colour="black",size=15),
        axis.text.x = element_text(colour="grey10",size=15),
        axis.title.x = element_text(colour="black",size=15),
        title = element_text(size=15),
        legend.text=element_text(size=15))+
  ggtitle("Database: Herring")


all_plot<-gridExtra::grid.arrange(masked_herring, masked_fish, nrow=2)
ggsave(plot=all_plot, filename="masked_proportion_by_region_herring_actinopterygii.pdf")

herrign_only<-table %>% 
  filter(Genotype!="NS") %>%
  filter(Database=="HerringCombinedDB") %>%
  ggplot(aes(y=Masked_Proportion, x=Chromosome, fill=Genotype))+
  geom_boxplot()+
  scale_fill_manual(values=c("#D95F02","#1F78B4")) +
  ylim(0, 45)+
  theme_classic()+
  ylab("Masked Proportion (%)")+
  theme(axis.text.y = element_text(colour="grey10",size=7),
        axis.title.y = element_text(colour="black",size=7),
        axis.text.x = element_text(colour="grey10",size=7),
        axis.title.x = element_text(colour="black",size=7),
        title = element_text(size=7),
        legend.text=element_text(size=7))
  ggtitle("Database: Herring")

ggsave(plot=herrign_only, filename="masked_proportion_by_region_herring_only.pdf", height=6, width=8)
ggsave(plot=herrign_only, filename="masked_proportion_by_region_herring_only_2023-04-21.pdf", height=40, width=200, units="mm")
ggsave(plot=herrign_only, filename="masked_proportion_by_region_herring_only_2023-04-02.pdf", height=40, width=200, units="mm")

## Test ####
## Chr 6 ####
chr6_table<-table %>% filter(Database=="HerringCombinedDB" & Chromosome=="Chr6")

games_howell_test(chr6_table, Masked_Proportion ~ Genotype, conf.level = 0.95, detailed =T)
# .y.               group1 group2    n1    n2 estimate conf.low conf.high    se statistic    df p.adj p.adj.signif method      
# * <chr>             <chr>  <chr>  <int> <int>    <dbl>    <dbl>     <dbl> <dbl>     <dbl> <dbl> <dbl> <chr>        <chr>       
#   1 Masked_Proportion S      N          6     6  -0.0700   -0.616     0.476 0.172     0.289  9.35 0.779 ns           Games-Howell

## Chr 12 ####
chr12_table<-table %>% filter(Database=="HerringCombinedDB" & Chromosome=="Chr12")

games_howell_test(chr12_table, Masked_Proportion ~ Genotype, conf.level = 0.95, detailed =T)
# A tibble: 1 × 14
# .y.               group1 group2    n1    n2 estimate conf.low conf.high    se statistic    df p.adj p.adj.signif method      
# * <chr>             <chr>  <chr>  <int> <int>    <dbl>    <dbl>     <dbl> <dbl>     <dbl> <dbl> <dbl> <chr>        <chr>       
#   1 Masked_Proportion S      N          6     6     1.25    -2.05      4.55  1.02     0.864  8.56 0.411 ns           Games-Howell


## Chr 17 ####
chr17_table<-table %>% filter(Database=="HerringCombinedDB" & Chromosome=="Chr17") %>% filter(Genotype!="NS")

games_howell_test(chr17_table, Masked_Proportion ~ Genotype, conf.level = 0.95, detailed =T)
# # A tibble: 1 × 14
# .y.               group1 group2    n1    n2 estimate conf.low conf.high    se statistic    df p.adj p.adj.signif method      
# * <chr>             <chr>  <chr>  <int> <int>    <dbl>    <dbl>     <dbl> <dbl>     <dbl> <dbl> <dbl> <chr>        <chr>       
#   1 Masked_Proportion S      N          5     6    0.119    -1.87      2.11 0.572     0.147  5.85 0.888 ns           Games-Howell

## Chr 23 ####
chr23_table<-table %>% filter(Database=="HerringCombinedDB" & Chromosome=="Chr23") %>% filter(Genotype!="NS")

games_howell_test(chr23_table, Masked_Proportion ~ Genotype, conf.level = 0.95, detailed =T)
# # A tibble: 1 × 14
# .y.               group1 group2    n1    n2 estimate conf.low conf.high    se statistic    df p.adj p.adj.signif method      
# * <chr>             <chr>  <chr>  <int> <int>    <dbl>    <dbl>     <dbl> <dbl>     <dbl> <dbl> <dbl> <chr>        <chr>       
#   1 Masked_Proportion S      N          5     6     5.41    0.223      10.6  1.62      2.36  9.00 0.043 *            Games-Howell

table_all<-table %>% filter(Genotype!="NS" & Database=="HerringCombinedDB")
table_all$group<-paste0(table_all$Genotype, "_", table_all$Chromosome)
table_all$group<-factor(table_all$group,levels = c("N_Genome", "S_Chr6", "N_Chr6", "S_Chr12", "N_Chr12", "S_Chr17", "N_Chr17", "S_Chr23", "N_Chr23"))

all_results<-data.frame(games_howell_test(table_all, Masked_Proportion ~ group, conf.level = 0.95, detailed =F))
write.table(all_results[c(1,2,3,4,5,6,7,8,9,22,31,36),], file="games_howell_test_results.txt", col.names=T, row.names = F, quote = F, sep = "\t")

all_results_wilcox<-data.frame(wilcox_test(table_all, Masked_Proportion ~ group))
all_results_wilcox[c(1,2,3,4,5,6,7,8,9,22,31,36),]


# select the results we want (against genome, and within inversion)
all_results_target<-all_results[c(1,2,3,4,5,6,7,8,9,36,31,22),]

annotation_df <- data.frame()

# code from here: https://www.r-bloggers.com/2017/06/add-p-values-and-significance-levels-to-ggplots/
# https://github.com/const-ae/ggsignif

my_comparisons <- list( c("N_Chr12", "S_Chr12"), c("N_Chr12", "Genome"), c("Genome", "S_Chr12"),
                        c("N_Chr6", "S_Chr6"), c("N_Chr6", "Genome"), c("Genome", "S_Chr6"),
                        c("N_Chr17", "S_Chr17"), c("N_Chr17", "Genome"), c("Genome", "S_Chr17"),
                        c("N_Chr23", "S_Chr23"), c("N_Chr23", "Genome"), c("Genome", "S_Chr23"))

test_result<-ggboxplot(table_all, x = "group", y = "Masked_Proportion",
           palette = "jco")+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 90)  

ggsave(plot=test_result, filename="boxplot_with_test_results.txt")

