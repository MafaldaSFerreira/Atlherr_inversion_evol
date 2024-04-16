#install.packages("data.table")
#install.packages("devtools")
#devtools::install_github("petrelharp/local_pca/lostruct")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#install.packages("optparse")
#BiocManager::install("Rsamtools")
#library(devtools)
#install_github("petrelharp/templater")

####################################
library(lostruct)
library(data.table)
library(ggplot2)
library(tidyverse)
####################################

####################################
my_thm= list( theme_bw(),
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme(panel.border = element_blank(),legend.position = "none"),
              theme(axis.line = element_line(colour = "black")),
              theme(axis.text.y = element_text(colour="grey10",size=7),
                    axis.title.y = element_text(colour="black",size=7),
                    axis.text.x = element_text(colour="grey10",size=7),
                    axis.title.x = element_text(colour="black",size=7)))
####################################

# Set the PATH otherwise read_vcf will not find bcftools!
#export_PATH=/opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin
#Sys.setenv(PATH="/opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin")

# Set working directory:
#setwd("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/LocalPCA/F2_mis20_maf0.01_29032023/dataset_35_individuals")
#save.image("/Users/maffe217/Dropbox/Mac/Documents/Postdoc/Repositories/Atlherr_inversion_evol_RData/4ii.LocalPCA/localPCA.RData")

#chr6
mds_results<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_0008/mds_coords.csv")
coords<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_0008/chr6.regions.csv")
pcas<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_0008/chr6.pca.csv")

cbind(pcas,coords)->results
results$mid<-(results$end+results$start)/2

# chr6
s=22282765
e=24868581
chr="chr6"

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

### Chr 6 ####
p1<-results %>%
  ggplot() +
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=-0.7,ymax=0.7),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  
  geom_line(aes(x=mid/1000000,y=PC_1_F1),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F2),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F3),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F4),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F5),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F6),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn3_Fehmarn_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn44_Fehmarn_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn6_Fehmarn_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle100_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle54_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle98_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF16_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF18_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF19_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF21_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM14_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM15_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM16_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM19_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  
  geom_line(aes(x=mid/1000000,y=PC_1_AAL1_CelticSea_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AAL2_CelticSea_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AAL3_Celticsea_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS10),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS4),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS5),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS7),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS8),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS2),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK1_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK2_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK3_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z12_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z14_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z4_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  
  ylim(-0.7,0.7)+
  xlim(s/1000000-2,e/1000000+2)+
  xlab("Position (Mb)") +
  ylab("PC1") +
  my_thm

ggsave(plot=p1,height=1.3,width=4,dpi=300, filename="figures/chr6_all_SUK_all_Baltic.PC1_1.3x4.pdf", useDingbats=FALSE)

### Chr 12 ######
# chr 12
s=17826318
e=25603093
chr="chr12"

mds_results<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_0007/mds_coords.csv")
coords<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_0007/chr12.regions.csv")
pcas<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_0007/chr12.pca.csv")

cbind(pcas,coords)->results
results$mid<-(results$end+results$start)/2

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

# chromosome 12
p1<-results %>%
  ggplot() +
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=-0.7,ymax=0.7),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(aes(x=mid/1000000,y=PC_1_AAL1_CelticSea_Atlantic_Winter),color='black')+

  geom_line(aes(x=mid/1000000,y=PC_1_F1),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F2),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F3),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F4),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F5),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F6),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn3_Fehmarn_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn44_Fehmarn_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn6_Fehmarn_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle100_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle54_Gavle_Baltic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle98_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF16_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF18_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF19_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF21_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM14_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM15_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM16_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM19_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  
  geom_line(aes(x=mid/1000000,y=PC_1_AAL2_CelticSea_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AAL3_Celticsea_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS10),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS4),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS5),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS7),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS8),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS2),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK1_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK2_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK3_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z12_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z14_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z4_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  
  ylim(-0.7,0.7)+
  xlim(s/1000000-2,e/1000000+2)+
  xlab("Position (Mb)") +
  ylab("PC1") +
  my_thm

#ggsave(plot=p1,height=3,width=10,dpi=300, filename="figures/chr12_all_SUK_all_Baltic.PC1.pdf", useDingbats=FALSE)
#ggsave(plot=p1,height=2,width=5,dpi=300, filename="figures/chr12_all_SUK_all_Baltic.PC1.2x5.pdf", useDingbats=FALSE)
ggsave(plot=p1,height=1.3,width=4,dpi=300, filename="figures/chr12_all_SUK_all_Baltic.PC1.1.3x4.pdf", useDingbats=FALSE)

### Chr 17 ######

mds_results<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_0009/mds_coords.csv")
coords<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_0009/chr17.regions.csv")
pcas<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_0009/chr17.pca.csv")

# chr17
s=25805445
e=27568511
chr="chr17"

cbind(pcas,coords)->results
results$mid<-(results$end+results$start)/2

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))
# chromosome 17
p1<-results %>%
  ggplot() +
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=-0.7,ymax=0.7),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
 
  geom_line(aes(x=mid/1000000,y=PC_1_F5),color='black',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F6),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn3_Fehmarn_Baltic_Autumn),color='black',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn44_Fehmarn_Baltic_Autumn),color='black',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn6_Fehmarn_Baltic_Autumn),color='black',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F1),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F2),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F3),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F4),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle100_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle54_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle98_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF16_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF18_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF19_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF21_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM14_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM15_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM16_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM19_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  
  geom_line(aes(x=mid/1000000,y=PC_1_AAL1_CelticSea_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AAL2_CelticSea_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AAL3_Celticsea_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS10),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS4),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS5),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS7),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS8),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS2),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK1_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK2_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK3_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z12_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z14_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z4_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  
  ylim(-0.7,0.7)+
  xlim(s/1000000-2,e/1000000)+
  xlab("Position (Mb)") +
  ylab("PC1") +
  my_thm

ggsave(plot=p1,height=3,width=10,dpi=300, filename="figures/chr17_all_SUK_all_Baltic.PC1.pdf", useDingbats=FALSE)
ggsave(plot=p1,height=2,width=5,dpi=300, filename="figures/chr17_all_SUK_all_Baltic.PC1.2x5.pdf", useDingbats=FALSE)
ggsave(plot=p1,height=1.3,width=4,dpi=300, filename="figures/chr17_all_SUK_all_Baltic.PC1.1.3x4.pdf", useDingbats=FALSE)

## CHR23 #####

s=16226443
e=17604273
chr="chr23"

mds_results<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_00010/mds_coords.csv")
coords<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_00010/chr23.regions.csv")
pcas<-read.csv("lostruct_results/type_snp_size_200_weights_none_jobid_00010/chr23.pca.csv")


cbind(pcas,coords)->results
results$mid<-(results$end+results$start)/2

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

p1<-results %>%
  ggplot() +
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=-0.7,ymax=0.7),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  
 
  geom_line(aes(x=mid/1000000,y=PC_1_F5),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F6),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn3_Fehmarn_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn44_Fehmarn_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Fehmarn6_Fehmarn_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F1),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  
  
  geom_line(aes(x=mid/1000000,y=PC_1_F2),color='black',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F3),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_F4),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle100_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle54_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Gavle98_Gavle_Baltic_Autumn),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF16_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF18_HastKar_Baltic_Spring),color='black',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM19_HastKar_Baltic_Spring),color='black',alpha=0.5, linewidth=0.5)+
  
  geom_line(aes(x=mid/1000000,y=PC_1_BF21_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM14_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM15_HastKar_Baltic_Spring),color='black',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BM16_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_BF19_HastKar_Baltic_Spring),color='#1F78B4',alpha=0.5, linewidth=0.5)+
  
  geom_line(aes(x=mid/1000000,y=PC_1_AAL1_CelticSea_Atlantic_Winter),color='black',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AAL2_CelticSea_Atlantic_Winter),color='black',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AAL3_Celticsea_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS10),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS4),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS5),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS7),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS8),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_CS2),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK1_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK2_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_AK3_Downs_Atlantic_Winter),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z12_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z14_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  geom_line(aes(x=mid/1000000,y=PC_1_Z4_IsleofMan_Atlantic_Autumn),color='#D95F02',alpha=0.5, linewidth=0.5)+
  
  ylim(-0.7,0.7)+
  xlim(s/1000000-2,e/1000000+2)+
  xlab("Position (Mb)") +
  ylab("PC1") +
  my_thm

ggsave(plot=p1,height=1.3,width=4,dpi=300, filename="figures/chr23_all_SUK_all_Baltic.PC1.1.3x4.pdf", useDingbats=FALSE)


#### SessionInfo ####
# R version 4.3.0 (2023-04-21)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS 14.4.1

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# time zone: Europe/Stockholm
# tzcode source: internal

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] data.table_1.14.8   lostruct_0.0.0.9000 ggrastr_1.0.2       cowplot_1.1.1       lubridate_1.9.3    
#  [6] forcats_1.0.0       stringr_1.5.0       dplyr_1.1.3         purrr_1.0.2         readr_2.1.4        
# [11] tidyr_1.3.0         tibble_3.2.1        ggplot2_3.4.4       tidyverse_2.0.0    

# loaded via a namespace (and not attached):
#  [1] gtable_0.3.4      compiler_4.3.0    tidyselect_1.2.0  ggbeeswarm_0.7.2  gridExtra_2.3     scales_1.2.1     
#  [7] R6_2.5.1          labeling_0.4.3    generics_0.1.3    Cairo_1.6-1       munsell_0.5.0     pillar_1.9.0     
# [13] tzdb_0.4.0        rlang_1.1.1       utf8_1.2.3        stringi_1.7.12    timechange_0.2.0  cli_3.6.1        
# [19] withr_2.5.1       magrittr_2.0.3    grid_4.3.0        rstudioapi_0.15.0 hms_1.1.3         beeswarm_0.4.0   
# [25] lifecycle_1.0.3   vipor_0.4.5       vctrs_0.6.4       glue_1.6.2        farver_2.1.1      fansi_1.0.5      
# [31] colorspace_2.1-0  tools_4.3.0       pkgconfig_2.0.3  