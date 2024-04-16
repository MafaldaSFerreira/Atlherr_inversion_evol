# Recombination heatmaps
# Mafalda S. Ferreira, UU, Uppsala

# Packages ####
library(ggrastr)
library(viridis)
library(cowplot)
library(gridExtra)
library(tidyverse)

# Working directory ####
setwd("~/Documents/Postdoc/Project_Herring/Inversion_Project/Recombination/repetition_miss20_maf0.1_03032023/datafiles/R2_vcf.snps.miss90GQ30")
save.image("~/Documents/Postdoc/Repositories/Atlherr_inversion_evol_RData/4iv.recombination/recombination_heatmaps.RData")

# Supplementary Figure 11 - Combined R2 heatmaps #####

#### Chr6 ####
##### All ####
homozygotes_chr6 <- read.table("only_homozygotes_chr6.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.ld",header=T)

homozygotes_chr6<-na.omit(homozygotes_chr6)

homozygotes_chr6_plot<-ggplot(homozygotes_chr6) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  #scale_fill_continuous(low="white", high="black")+
  scale_fill_viridis(option = "B", direction=-1)+
  #ylim(21.3,25)+
  #xlim(21.3,25)+
  xlab("Chr6 (Mb)")+
  ylab("NN+SS")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### South ####
south_chr6 <- read.table("homozygotes_South_chr6.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2.geno.ld",header=T)

south_chr6<-na.omit(south_chr6)

south_chr6_plot<-ggplot(south_chr6) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  #scale_fill_continuous(low="white", high="black")+
  scale_fill_viridis(option = "B", direction=-1)+
  #ylim(21.3,25)+
  #xlim(21.3,25)+
  xlab("Chr6 (Mb)")+
  ylab("SS")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### North ####
north_chr6 <- read.table("homozygotes_North_chr6.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2.geno.ld",header=T)

north_chr6<-na.omit(north_chr6)

north_chr6_plot<-ggplot(north_chr6) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  #scale_fill_continuous(low="white", high="black")+
  scale_fill_viridis(option = "B", direction=-1)+
  #ylim(21.3,25)+
  #xlim(21.3,25)+
  xlab("Chr6 (Mb)")+
  ylab("NN")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

#### Plot All Chr6 ####
all<-gridExtra::grid.arrange(homozygotes_chr6_plot, north_chr6_plot, south_chr6_plot, nrow=3)
ggsave(all, filename ="figures/chr6.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.pdf", dpi=300, unit="mm", height=110, width=160)

#### Chr12 ####
##### All ####
homozygotes_chr12 <- read.table("only_homozygotes_chr12.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.ld",header=T)

homozygotes_chr12<-na.omit(homozygotes_chr12)

homozygotes_chr12_plot<-ggplot(homozygotes_chr12) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  scale_fill_viridis(option = "B", direction=-1)+
  xlab("Chr12 (Mb)")+
  ylab("NN + SS")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### South ####
south_chr12 <- read.table("homozygotes_South_chr12.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2.geno.ld",header=T)

south_chr12<-na.omit(south_chr12)

south_chr12_plot<-ggplot(south_chr12) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  scale_fill_viridis(option = "B", direction=-1)+
  xlab("Chr12 (Mb)")+
  ylab("SS")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### North ####
north_chr12 <- read.table("homozygotes_North_chr12.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2.geno.ld",header=T)

north_chr12<-na.omit(north_chr12)

north_chr12_plot<-ggplot(north_chr12) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  #scale_fill_continuous(low="white", high="black")+
  scale_fill_viridis(option = "B", direction=-1)+
  #ylim(21.3,25)+
  #xlim(21.3,25)+
  xlab("Chr12 (Mb)")+
  ylab("NN")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### Plot All Chr 12####
all_chr12<-gridExtra::grid.arrange(homozygotes_chr12_plot, north_chr12_plot, south_chr12_plot, nrow=3)

ggsave(all_chr12, filename ="figures/chr12.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.png", dpi=300, height=6, width=6)
ggsave(all_chr12, filename ="figures/chr12.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.pdf", dpi=300, units="mm", height=110, width=160)

#### Chr17 ####
##### ALL ####
homozygotes_chr17 <- read.table("only_homozygotes_chr17.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.ld",header=T)

homozygotes_chr17<-na.omit(homozygotes_chr17)

homozygotes_chr17_plot<-ggplot(homozygotes_chr17) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  scale_fill_viridis(option = "B", direction=-1)+
  xlab("Chr17 (Mb)")+
  ylab("SS + NN")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### South ####
south_chr17 <- read.table("homozygotes_South_chr17.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2.geno.ld",header=T)

south_chr17<-na.omit(south_chr17)

south_chr17_plot<-ggplot(south_chr17) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  scale_fill_viridis(option = "B", direction=-1)+
  xlab("Chr17 (Mb)")+
  ylab("SS")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### North ####
north_chr17 <- read.table("homozygotes_North_chr17.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2.geno.ld",header=T)

north_chr17<-na.omit(north_chr17)

north_chr17_plot<-ggplot(north_chr17) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)),dpi=300)+
  scale_fill_viridis(option = "B", direction=-1)+
  xlab("Chr17 (Mb)")+
  ylab("NN")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### Plot All Chr17 ####
all_chr17<-gridExtra::grid.arrange(homozygotes_chr17_plot, north_chr17_plot, south_chr17_plot, nrow=3)

ggsave(all_chr17, filename ="figures/chr17.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.png", dpi=300, height=6, width=6)
ggsave(all_chr17, filename ="figures/chr17.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.pdf", dpi=300, units="mm", height=110, width=160)

#### Chr23 ####
##### ALL ####
homozygotes_chr23 <- read.table("only_homozygotes_chr23.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.ld",header=T)

homozygotes_chr23<-na.omit(homozygotes_chr23)

homozygotes_chr23_plot<-ggplot(homozygotes_chr23) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  scale_fill_viridis(option = "B", direction=-1)+
  xlab("Chr23 (Mb)")+
  ylab("SS + NN")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### South ####
south_chr23 <- read.table("homozygotes_South_chr23.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2.geno.ld",header=T)

south_chr23<-na.omit(south_chr23)

south_chr23_plot<-ggplot(south_chr23) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  scale_fill_viridis(option = "B", direction=-1)+
  xlab("Chr23 (Mb)")+
  ylab("SS")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### North ####
north_chr23 <- read.table("homozygotes_North_chr23.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2.geno.ld",header=T)

north_chr23<-na.omit(north_chr23)

north_chr23_plot<-ggplot(north_chr23) +
  rasterise(geom_tile(aes(x=POS1/1000000, y=POS2/1000000, fill=R.2)), dpi=300)+
  scale_fill_viridis(option = "B", direction=-1)+
  xlab("Chr23 (Mb)")+
  ylab("NN")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10))

##### Plot All Chr23 ####
all_chr23<-gridExtra::grid.arrange(homozygotes_chr23_plot, north_chr23_plot, south_chr23_plot, nrow=3)

ggsave(all_chr23, filename ="figures/chr23.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.png", dpi=300, height=6, width=6)
ggsave(all_chr23, filename ="figures/chr23.snps.miss90.GQ30.thin5kb.noLdwd.R2.geno.pdf", dpi=300, units="mm", height=110, width=160)

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
