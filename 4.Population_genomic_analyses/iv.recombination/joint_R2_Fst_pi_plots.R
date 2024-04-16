# Join Fst, Pi and R2 plots for Figure 5
# Mafalda S. Ferreira, UU, Uppsala

# Packages ####
library(ggrastr)
library(viridis)
library(cowplot)
library(gridExtra)
library(tidyverse)

# Working directory ####
setwd("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.1_03032023")

# Figure 5 Plots A-D #####

#### Chr6 ####
s=22282765
e=24868581
chr="chr6"
POP1=c("CS_chr6")
POP2=c("BS_chr6")
XMinLim=21.3
XMaxLim=26

table_fst<-read.table(file="~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/datafiles/chr6_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_fst.txt", header=T)
table_pi<-read.table(file="~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/datafiles/chr6_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)

LD.raw <- read.table("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Recombination/repetition_miss20_maf0.1_03032023/datafiles/R2_vcf.snps.miss90GQ30_inversion/only_homozygotes_chr6.snps.miss90.GQ30.thin5kb.noLdwd.inversion.R2.geno.ld",header=T)

#### Chr12 ####

# chr 12
s=17826318
e=25603093
chr="chr12"
POP1=c("CS_chr12")
POP2=c("BS_chr12")
XMinLim=s/1000000-1
XMaxLim=e/1000000+1

table_fst<-read.table(file="~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/datafiles/chr12_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_fst.txt", header=T)
table_pi<-read.table(file="~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/datafiles/chr12_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)
LD.raw <- read.table("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Recombination/repetition_miss20_maf0.1_03032023/datafiles/R2_vcf.snps.miss90GQ30_inversion/only_homozygotes_chr12.snps.miss90.GQ30.thin5kb.noLdwd.inversion.R2.geno.ld",header=T)

#### Chr17 ####

s=25805445
e=27568511
chr="chr17"
POP1=c("CS_chr17")
POP2=c("BS_chr17")
XMinLim=s/1000000-1
XMaxLim=e/1000000

table_fst<-read.table(file="~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/datafiles/chr17_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_fst.txt", header=T)
table_pi<-read.table(file="~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/datafiles/chr17_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)
LD.raw <- read.table("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Recombination/repetition_miss20_maf0.1_03032023/datafiles/R2_vcf.snps.miss90GQ30_inversion/only_homozygotes_chr17.snps.miss90.GQ30.thin5kb.noLdwd.inversion.R2.geno.ld",header=T)


#### Chr23 ####
s=16226443
e=17604273
chr="chr23"
POP1=c("CS_chr23")
POP2=c("F_chr23")
XMinLim=s/1000000-1
XMaxLim=e/1000000+1

table_fst<-read.table(file="~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/datafiles/chr23_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_fst.txt", header=T)
table_pi<-read.table(file="~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/datafiles/chr23_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)
LD.raw <- read.table("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Recombination/repetition_miss20_maf0.1_03032023/datafiles/R2_vcf.snps.miss90GQ30_inversion/only_homozygotes_chr23.snps.miss90.GQ30.thin5kb.noLdwd.inversion.R2.geno.ld",header=T)


#### Plots ####
table_fst$mid<-(table_fst$window_pos_1 + table_fst$window_pos_2 )/2
fst_pixy<-table_fst %>%
  filter(chromosome== chr & window_pos_1 >= s-1000000 & window_pos_2 <=e+1000000) %>%
  filter(pop1==POP1 & pop2==POP2)

table_pi$mid<-(table_pi$window_pos_1 + table_pi$window_pos_2 )/2
pi_pixy<-table_pi %>%
  filter(no_sites>=10000) %>%
  filter(chromosome== chr & window_pos_1 >= s-1000000 & window_pos_2 <=e+1000000) %>%
  filter(pop==POP1 | pop==POP2)

LD.raw<-na.omit(LD.raw)

p1<-ggplot(LD.raw) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  #scale_fill_gradient2(low="white", high="black")+
  scale_fill_viridis(option = "A", direction=-1)+
  ylim(XMinLim,XMaxLim)+
  xlim(XMinLim,XMaxLim)+
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " position (Mb)"))+
  ylab(paste0("Chr ", str_remove(chr, "chr"), " (Mb)"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(), legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_text(colour="black",size=7),
        axis.text.x = element_text(colour="black",size=7),
        axis.title.x = element_text(colour="black",size=7))


candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

quant_0.99_fst<-table_fst %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  summarise(quantile_CS_F=quantile(avg_wc_fst, 0.99, na.rm=T))

p2<-ggplot(fst_pixy)+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=1),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(aes(x=mid/1000000,y=avg_wc_fst),color="black", linewidth=1)+
  geom_hline(yintercept = quant_0.99_fst$quantile_CS_F, color="black",linetype="dashed")+
  #xlim(s/1000000-1,e/1000000+1)+
  xlim(XMinLim,XMaxLim)+
  ylab("Fst") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="grey10",size=7),
        axis.title.y = element_text(colour="black",size=7),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

quant_0.01_pi<-table_pi %>%
  filter(pop!="Vancouver") %>%
  summarise(quantile_CS_F=quantile(avg_pi, 0.05, na.rm=T))


p3<-pi_pixy %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.01),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(aes(x=mid/1000000, y=avg_pi, color=pop), size=1, alpha=0.7)+
  scale_color_manual(values=c("#1F78B4", "#D95F02"))+
  geom_hline(yintercept = quant_0.01_pi$quantile_CS_F, color="black",linetype="dashed")+
  xlim(XMinLim,XMaxLim)+
  ylim(0,0.01)+
  ylab("Pi")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="grey10",size=7),
        axis.title.y = element_text(colour="black",size=7),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())


#### Save Plots ####
all_PLOTS<-plot_grid(p2, p3, p1, align = "v", nrow = 3)

ggsave(all_PLOTS, filename=paste0("figures/Fst_PI_R2_",chr,"_4X3_withscale.pdf"), dpi=300, width=4, height=3)

# Session Info #
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
#   [1] lubridate_1.9.3   forcats_1.0.0     stringr_1.5.0     dplyr_1.1.3       purrr_1.0.2       readr_2.1.4      
# [7] tidyr_1.3.0       tibble_3.2.1      ggplot2_3.4.4     tidyverse_2.0.0   gridExtra_2.3     cowplot_1.1.1    
# [13] viridis_0.6.4     viridisLite_0.4.2 ggrastr_1.0.2    
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.4        crayon_1.5.2        compiler_4.3.0      BiocManager_1.30.22 tidyselect_1.2.0   
# [6] ggbeeswarm_0.7.2    scales_1.2.1        R6_2.5.1            labeling_0.4.3      generics_0.1.3     
# [11] Cairo_1.6-1         munsell_0.5.0       pillar_1.9.0        tzdb_0.4.0          rlang_1.1.1        
# [16] utf8_1.2.3          stringi_1.7.12      timechange_0.2.0    cli_3.6.1           withr_2.5.1        
# [21] magrittr_2.0.3      grid_4.3.0          rstudioapi_0.15.0   hms_1.1.3           beeswarm_0.4.0     
# [26] lifecycle_1.0.3     vipor_0.4.5         vctrs_0.6.4         glue_1.6.2          farver_2.1.1       
# [31] fansi_1.0.5         colorspace_2.1-0    tools_4.3.0         pkgconfig_2.0.3    
