####################################################
# Genotype heatmaps
# Mafalda Ferreira,
# Uppsala, Sweden
####################################################
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrastr)

session<-sessionInfo()

###
# R version 4.3.0 (2023-04-21)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS 14.4.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Europe/Stockholm
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrastr_1.0.2   cowplot_1.1.1   lubridate_1.9.3 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.3     purrr_1.0.2    
# [8] readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.4      compiler_4.3.0    tidyselect_1.2.0  ggbeeswarm_0.7.2  gridExtra_2.3     scales_1.2.1     
# [7] R6_2.5.1          generics_0.1.3    munsell_0.5.0     pillar_1.9.0      tzdb_0.4.0        rlang_1.1.1      
# [13] utf8_1.2.3        stringi_1.7.12    timechange_0.2.0  cli_3.6.1         withr_2.5.1       magrittr_2.0.3   
# [19] grid_4.3.0        rstudioapi_0.15.0 hms_1.1.3         beeswarm_0.4.0    lifecycle_1.0.3   vipor_0.4.5      
# [25] vctrs_0.6.4       glue_1.6.2        fansi_1.0.5       colorspace_2.1-0  tools_4.3.0       pkgconfig_2.0.3  
###


#setwd("~/Documents/Postdoc/Project_Herring/Inversion_Project/Genotypes/F2_miss20maf0.1_27022023/")
#save.image("~/Documents/Postdoc/Repositories/Atlherr_inversion_evol_RData/4i.genotype_heatmap/genotype_heatmap.RData")

table12<-read.table("herring_sentieon_91ind_190521.chr12.NSvNE_logp15.nomono.nomulti.sep.genotypes")
table6<-read.table("herring_sentieon_91ind_190521.chr6.NSvNE_logp15.nomono.nomulti.sep.genotypes")
table17<-read.table("herring_sentieon_91ind_190521.chr17.NSvNE_logp15.nomono.nomulti.sep.genotypes")
table23<-read.table("herring_sentieon_91ind_190521.chr23.NSvNE_logp15.nomono.nomulti.sep.genotypes")


hd<-read.table("../header.txt")
hd[,c("V31", "V32", "V33", "V36", "V41", "V42")]<-c("BS1", "BS2", "BS3", "BS4", "BS5", "BS6")

colnames(table6)<-hd[1,]
colnames(table12)<-hd[1,]
colnames(table17)<-hd[1,]
colnames(table23)<-hd[1,]

indv<-read.table("../individuals.txt")

# For publication, we want to change the names of F individuls to BS...
indv[c(19,20,21,22,23,24),]<-c("BS1", "BS2", "BS3", "BS4", "BS5", "BS6")

table6[table6=="0/0"]<-0
table6[table6=="0/1"]<-1
table6[table6=="1/0"]<-1
table6[table6=="1/1"]<-2
table6[table6=="./."]<-NA

table12[table12=="0/0"]<-0
table12[table12=="0/1"]<-1
table12[table12=="1/0"]<-1
table12[table12=="1/1"]<-2
table12[table12=="./."]<-NA

table17[table17=="0/0"]<-0
table17[table17=="0/1"]<-1
table17[table17=="1/0"]<-1
table17[table17=="1/1"]<-2
table17[table17=="./."]<-NA

table23[table23=="0/0"]<-0
table23[table23=="0/1"]<-1
table23[table23=="1/0"]<-1
table23[table23=="1/1"]<-2
table23[table23=="./."]<-NA


# solution
table6 %>% 
  drop_na %>% # forgot to remove missing data
  pivot_longer(cols=3:93,names_to="individuals",values_to="genotypes") -> table_genos_6

table12 %>% 
  drop_na %>% # forgot to remove missing data
  pivot_longer(cols=3:93,names_to="individuals",values_to="genotypes") -> table_genos_12

table17 %>% 
  drop_na %>% # forgot to remove missing data
  pivot_longer(cols=3:93,names_to="individuals",values_to="genotypes") -> table_genos_17

table23 %>% 
  drop_na %>% # forgot to remove missing data
  pivot_longer(cols=3:93,names_to="individuals",values_to="genotypes") -> table_genos_23

# polarize chr 17 and chr23:
# We need to polarize the genotypes in chr23 and chr17.
# let's use the bottom individuals as reference:
reference_indvs<-c("CS2", "CS8", "CS7", "CS5", "CS4", "CS10",
                   "AAL1_CelticSea_Atlantic_Winter",
                   "AAL2_CelticSea_Atlantic_Winter",
                   "AAL3_Celticsea_Atlantic_Winter",
                   "AK1_Downs_Atlantic_Winter",
                   "AK2_Downs_Atlantic_Winter",
                   "AK3_Downs_Atlantic_Winter",
                   "Z12_IsleofMan_Atlantic_Autumn",
                   "Z14_IsleofMan_Atlantic_Autumn",
                   "Z4_IsleofMan_Atlantic_Autumn"
)

individuals<-names(table17_nomiss)[!names(table17_nomiss) %in% c("CHROM", "POS")]

# find rows where these individuals don't have a majority of "2" alleles.
table17_nomiss<-table17 %>% drop_na
table17_nomiss[,2:93]<-sapply(table17_nomiss[,c(2:93)], as.integer)

flip_rows <-rowSums(table17_nomiss[,reference_indvs])/(length(reference_indvs)*2) > 0.5

# Flip allele values for those rows
table17_nomiss[flip_rows, individuals]<-(table17_nomiss[flip_rows, individuals]-2)*-1

# this looks good. let's do it for chr23
table23_nomiss<-table23 %>% drop_na
table23_nomiss[,2:93]<-sapply(table23_nomiss[,c(2:93)], as.integer)

flip_rows <-rowSums(table23_nomiss[,reference_indvs])/(length(reference_indvs)*2) > 0.5

table23_nomiss[flip_rows, individuals]<-(table23_nomiss[flip_rows, individuals]-2)*-1

# Plots

table_genos_6$POS<-as.character(table_genos_6$POS)
table_genos_12$POS<-as.character(table_genos_12$POS)
table17_nomiss$POS<-as.character(table17_nomiss$POS)
table23_nomiss$POS<-as.character(table23_nomiss$POS)


chr6<-table_genos_6 %>%
  #filter(CHROM=="chr6") %>%
  #filter(genotypes=="0" | genotypes=="1" | genotypes=="2") %>%
  ggplot() +
  rasterise(geom_tile(aes(x=POS,y=factor(individuals,level=indv[,1]),fill=genotypes)), dpi=300)+
  scale_fill_manual(values=c("#1F78B4", "#CCCCCC", "#D95F02"))+
  xlab("Individuals")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour="black",size=7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

chr12<-table_genos_12 %>%
  ggplot() +
  rasterise(geom_tile(aes(x=POS,y=factor(individuals,level=indv[,1]),fill=genotypes)), dpi=300)+
  scale_fill_manual(values=c("#D95F02", "#CCCCCC", "#1F78B4"))+
  xlab("Individuals")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

chr17<-table17_nomiss %>%
  pivot_longer(cols=3:93,names_to="individuals",values_to="genotypes") %>%
  ggplot() +
  rasterise(geom_tile(aes(x=as.character(POS),y=factor(individuals,level=indv[,1]),fill=as.character(genotypes))), dpi=300)+
  scale_fill_manual(values=c("#D95F02", "#CCCCCC", "#1F78B4"))+
  xlab("Individuals")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

chr23<-table23_nomiss %>%
  pivot_longer(cols=3:93,names_to="individuals",values_to="genotypes") %>%
  ggplot() +
  rasterise(geom_tile(aes(x=as.character(POS),y=factor(individuals,level=indv[,1]),fill=as.character(genotypes))), dpi=300)+
  scale_fill_manual(values=c("#D95F02", "#CCCCCC", "#1F78B4"))+
  xlab("Individuals")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot_all<-plot_grid(chr6, chr12, chr17, chr23, align = "h", nrow = 1, rel_widths = c(0.35,0.5,0.1,0.05))

ggsave(plot = plot_all, filename = "figures/all_genotypes_combined_24-04-11.png", width = 8, height=8)
ggsave(plot = plot_all, filename = "figures/all_genotypes_combined_24-04-11.pdf", width = 8, height=8)
