### Plot summary statistics ###

### Libraries ####
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(biomaRt)
library(cowplot)

session<-sessionInfo()

save.image("~/Dropbox/Mac/Documents/Postdoc/Repositories/Atlherr_inversion_evol_RData/4iii.summary_stats/summary_stats.RData") # nolint: line_length_linter.

# Plot parameters ####
my_thm_bottom_plot= list( theme_bw(),
                          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
                          theme(panel.border = element_blank(),legend.position = "none"),
                          theme(axis.line = element_line(colour = "black")),
                          theme(axis.text.y = element_text(colour="black",size=7),
                                axis.title.y = element_text(colour="black",size=7),
                                axis.text.x = element_text(colour="black",size=7),
                                axis.title.x = element_text(colour="black",size=7)))

my_thm_middle_plots= list( theme_bw(),
                           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
                           theme(panel.border = element_blank(),legend.position = "none"),
                           theme(axis.line = element_line(colour = "black")),
                           theme(axis.text.y = element_text(colour="black",size=7),
                                 axis.title.y = element_text(colour="black",size=7),
                                 axis.text.x = element_blank(),
                                 axis.title.x = element_blank()))

my_thm3= list( theme_bw(),
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
               theme(panel.border = element_blank()),
               theme(axis.line = element_line(colour = "black")),
               theme(axis.text.y = element_text(colour="black",size=7),
                     axis.title.y = element_text(colour="black",size=7),
                     axis.text.x = element_text(colour="black",size=7),
                     axis.title.x = element_blank(),
                     legend.text=element_text(size=5),
                     legend.title = element_text(size=5)))



setwd("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023")

# Define plotting regions ####
#### Chr6 ####
s=22282765
e=24868581
chr="chr6"
POP1=c("CS_chr6")
POP2=c("BS_chr6")
XMinLim=s/1000000-1
XMaxLim=e/1000000+1

#Miss20Maf0.01
table_dxy<-read.table("datafiles/chr6_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")
# Miss20 Maf 0.01
table_fst<-read.table("datafiles/chr6_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_fst.txt", header=T)
# Miss 20 Maf 0.01
table_pi<-read.table("datafiles/chr6_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)

# create a subset table for the inversion
table_pi_inversion_chr6<-table_pi %>% filter(no_sites > 10000) %>% filter(pop==POP1 | pop==POP2) %>% 
  filter(chromosome==chr & window_pos_1 >=s & window_pos_2 <=e)

# create a subset table for the inversion
table_dxy_inversion_chr6<-table_dxy %>% filter(pop1==POP1 & pop2==POP2) %>% 
  filter(chromosome==chr & window_pos_1 >=s & window_pos_2 <=e)

write.table(table_dxy_inversion_chr6, file="datafiles/chr6_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.INVERSION.txt")

#### Chr12 ####

# chr 12
s=17826318
e=25603093
chr="chr12"
POP1=c("CS_chr12")
POP2=c("BS_chr12")
XMinLim=s/1000000-1
XMaxLim=e/1000000+1

#Miss20Maf0.01
table_dxy<-read.table("datafiles/chr12_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")
# Miss20 Maf 0.01
table_fst<-read.table("datafiles/chr12_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_fst.txt", header=T)
# Miss 20 Maf 0.01
table_pi<-read.table("datafiles/chr12_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)

table_pi_inversion_chr12<-table_pi %>% filter(no_sites > 10000) %>% filter(pop==POP1 | pop==POP2) %>% 
  filter(chromosome==chr & window_pos_1 >=s & window_pos_2 <=e)

write.table(table_pi_inversion_chr12, file="datafiles/chr12_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.INVERSION.txt")

# create a subset table for the inversion
table_dxy_inversion_chr12<-table_dxy %>% filter(pop1==POP1 & pop2==POP2) %>% 
  filter(chromosome==chr & window_pos_1 >=s & window_pos_2 <=e)


#### Chr17 ####

s=25805445
e=27568511
chr="chr17"
POP1=c("CS_chr17")
POP2=c("BS_chr17")
XMinLim=s/1000000-1
XMaxLim=e/1000000

#Miss20Maf0.01
table_dxy<-read.table("datafiles/chr17_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")
# Miss20 Maf 0.01
table_fst<-read.table("datafiles/chr17_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_fst.txt", header=T)
# Miss 20 Maf 0.01
table_pi<-read.table("datafiles/chr17_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)

table_pi_inversion_chr17<-table_pi %>% filter(no_sites > 10000) %>% filter(pop==POP1 | pop==POP2) %>% 
  filter(chromosome==chr & window_pos_1 >=s & window_pos_2 <=e)

write.table(table_pi_inversion_chr17, file="datafiles/chr17_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.INVERSION.txt")

# create a subset table for the inversion
table_dxy_inversion_chr17<-table_dxy %>% filter(pop1==POP1 & pop2==POP2) %>% 
  filter(chromosome==chr & window_pos_1 >=s & window_pos_2 <=e)


#### Chr23 ####
s=16226443
e=17604273
chr="chr23"
POP1=c("CS_chr23")
POP2=c("F_chr23")
XMinLim=s/1000000-1
XMaxLim=e/1000000+1

#Miss20Maf0.01
table_dxy<-read.table("datafiles/chr23_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")
# Miss20 Maf 0.01
table_fst<-read.table("datafiles/chr23_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_fst.txt", header=T)
# Miss 20 Maf 0.01
table_pi<-read.table("datafiles/chr23_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)

table_pi_inversion_chr23<-table_pi %>% filter(no_sites > 10000) %>% filter(pop==POP1 | pop==POP2) %>% 
  filter(chromosome==chr & window_pos_1 >=s & window_pos_2 <=e)

write.table(table_pi_inversion_chr23, file="datafiles/chr23_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.INVERSION.txt")

# create a subset table for the inversion
table_dxy_inversion_chr23<-table_dxy %>% filter(pop1==POP1 & pop2==POP2) %>% 
  filter(chromosome==chr & window_pos_1 >=s & window_pos_2 <=e)

# Main Figure 5 A-D ####

#### Divergence (Dxy) Not in Figure 5 ####
candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

table_dxy$mid<-(table_dxy$window_pos_1 + table_dxy$window_pos_2 )/2

quant_0.99<-table_dxy %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  summarise(quantile_CS_F=quantile(avg_dxy, 0.99, na.rm=T))

dxy_pixy<-table_dxy %>%
  filter(chromosome==chr) %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  filter(no_sites > 10000) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.02),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(aes(x=mid/1000000,y=avg_dxy),color="black",size=1)+
  geom_hline(yintercept = quant_0.99$quantile_CS_F,color="black",linetype="dashed")+
  xlim(XMinLim,XMaxLim)+
  ylim(0,0.02)+
  #ggtitle("pixy")+
  ylab("Dxy") +
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  my_thm_middle_plots

#### Differentiation Fst ####
table_fst$mid<-(table_fst$window_pos_1 + table_fst$window_pos_2 )/2

quant_0.99_fst<-table_fst %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  summarise(quantile_CS_F=quantile(avg_wc_fst, 0.99, na.rm=T))

fst_pixy<-table_fst %>%
  filter(chromosome== chr) %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  #filter(no_snps > 250) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=1),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(aes(x=mid/1000000,y=avg_wc_fst),color="black", size=1)+
  geom_hline(yintercept = quant_0.99_fst$quantile_CS_F, color="black",linetype="dashed")+
  xlim(XMinLim,XMaxLim)+
  #ggtitle("pixy")+
  ylab("Fst") +
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  my_thm_middle_plots

#### Nucleotide Diversity (Pi) ####
quant_0.99_pi<-table_pi %>%
  filter(pop!="Vancouver") %>%
  summarise(quantile_CS_F=quantile(avg_pi, 0.99, na.rm=T))

quant_0.01_pi<-table_pi %>%
  filter(pop!="Vancouver") %>%
  summarise(quantile_CS_F=quantile(avg_pi, 0.01, na.rm=T))

#mean_pi<-table_pi %>%
#  filter(pop!="Vancouver") %>%
#  summarise(quantile_CS_F=mean(avg_pi, na.rm=T))


table_pi$mid<-(table_pi$window_pos_1 + table_pi$window_pos_2 )/2

pi_pixy<-table_pi %>%
  filter(chromosome==chr) %>%
  filter(no_sites>=10000) %>%
  filter(pop!="Vancouver") %>%
  filter(pop==POP1 | pop==POP2) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.01),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(aes(x=mid/1000000, y=avg_pi, color=pop), size=1, alpha=0.7)+
  scale_color_manual(values=c("#D95F02", "#1F78B4"))+
  geom_hline(yintercept = quant_0.99_pi$quantile_CS_F, color="black",linetype="dashed")+
  xlim(XMinLim,XMaxLim)+
  ylim(0,0.01)+
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  ylab("Pi")+
  my_thm_middle_plots


#### Net Nucleotide Diversity (Da) ####
table_pi_10kb_no_sites <-table_pi %>% filter(no_sites>=10000)
table_dxy_10kb_no_sites <- table_dxy %>% filter(no_sites>=10000)

table_pi_CS<-table_pi_10kb_no_sites %>% filter(pop==POP1)
table_pi_F<-table_pi_10kb_no_sites %>% filter(pop==POP2)
table_dxy_CS_F<-table_dxy_10kb_no_sites %>% filter(pop1==POP1 & pop2==POP2)

table_pi_ALL<-left_join(table_pi_CS, table_pi_F, by=c("chromosome", "window_pos_1", "window_pos_2"),
                        suffix=c(".CS", ".BS"))

table_dxy_pi<-left_join(table_pi_ALL, table_dxy_CS_F, by=c("chromosome","window_pos_1", "window_pos_2"))

table_dxy_pi$da<- table_dxy_pi$avg_dxy - ((table_dxy_pi$avg_pi.CS + table_dxy_pi$avg_pi.BS)/2)

output_table<-table_dxy_pi %>%
  filter(chromosome==chr & window_pos_1>=s & window_pos_2<=e) %>%
  summarise( mean_da= mean(da*100, na.rm=T),
             median_da = median(da*100, na.rm=T), 
             mean_dxy = mean(avg_dxy*100, na.rm=T),
             avg_pi_CS = mean(avg_pi.CS*100, na.rm=T),
             avg_pi_BS = mean(avg_pi.BS*100, na.rm=T))

write.table(output_table, 
            file = paste0("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/summary_statistics/",chr,"_20kb_fst_pi_dxy_da_miss20maf0.01_pixy.tsv"),
            quote = F, row.names = F, col.names = T)

quant_0.99_da<-table_dxy_pi %>%
  summarise(quantile_CS_F=quantile(da, 0.99, na.rm=T))

da_pixy<-table_dxy_pi %>%
  filter(chromosome==chr) %>%
  filter(no_sites>=10000) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.005),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(aes(x=mid/1000000, y=da), size=1)+
  geom_hline(yintercept = quant_0.99_da$quantile_CS_F, color="black",linetype="dashed")+
  xlim(XMinLim,XMaxLim)+
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  ylab("da")+
  my_thm_bottom_plot


final_plot<- plot_grid(fst_pixy, pi_pixy, dxy_pixy, da_pixy, align = "v", nrow = 4)


ggsave(final_plot, 
       filename=paste0("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/figures/",chr,"_20kb_fst_pi_dxy_da_miss20maf0.01_pixy_4X6.pdf"), 
       dpi=300, width=4, height=6)


# Atlantic vs Pacific Herring ####

# Miss 20 maf 0.01
table_dxy<-read.table("datafiles/atlantic_vs_pacific_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")
table_pi<-read.table("datafiles/atlantic_vs_pacific_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt",header=T,sep="\t")

table_dxy %>%
  filter(no_sites>=10000) %>%
  summarise(mean_dxy=mean(avg_dxy, na.rm=T))


# miss 20 maf 0.01
mean_dxy
#1 0.004630083

table_pi_ATL<-subset(table_pi, table_pi$pop=="Atlantic")
table_pi_ATL$weigh_pi <- table_pi_ATL$avg_pi * table_pi_ATL$no_sites
sum(table_pi_ATL$weigh_pi, na.rm = T)/sum(table_pi_ATL$no_sites, na.rm=T)

table_pi_10kb_no_sites <-table_pi %>% filter(no_sites>=10000)
table_dxy_10kb_no_sites <- table_dxy %>% filter(no_sites>=10000)

table_pi_ATL<-table_pi_10kb_no_sites %>% filter(pop=="Atlantic")
table_pi_PAC<-table_pi_10kb_no_sites %>% filter(pop=="Pacific")
table_dxy_ALT_PAC<-table_dxy_10kb_no_sites %>% filter(pop1=="Atlantic" & pop2=="Pacific")

table_pi_ALL<-left_join(table_pi_ATL, table_pi_PAC, by=c("chromosome", "window_pos_1", "window_pos_2"),
                        suffix=c(".ATL", ".PAC"))

table_dxy_pi<-left_join(table_pi_ALL, table_dxy_ALT_PAC, by=c("chromosome","window_pos_1", "window_pos_2"))

table_dxy_pi$da<- table_dxy_pi$avg_dxy - ((table_dxy_pi$avg_pi.ATL + table_dxy_pi$avg_pi.PAC)/2)

output_table<-table_dxy_pi %>%
  #filter(chromosome=="chr1" | chromosome =="chr6" | chromosome =="chr12" | chromosome =="chr17" | chromosome =="chr23") %>%
  #filter(chromosome==chr & window_pos_1>=s & window_pos_2<=e) %>%
  summarise( mean_da= mean(da*100, na.rm=T),
             median_da = median(da*100, na.rm=T), 
             mean_dxy = mean(avg_dxy*100, na.rm=T),
             avg_pi_ATL = mean(avg_pi.ATL*100, na.rm=T),
             avg_pi_PAC = mean(avg_pi.PAC*100, na.rm=T))

write.table(output_table, 
            file = paste0("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/summary_statistics/Atlantic_vs_Herring_20kb_fst_pi_dxy_da_miss20maf0.01_pixy.tsv"),
            quote = F, row.names = F, col.names = T)

# Supplementary Fig 10A. Dxy between Northern and Southern homozygotes versus Atlantic vs Pacific herring  ####
# Generate inputs from MAF 0.01 and 
# Atlantic vs Pacific 

# Read input files:
#Miss20Maf0.01
table_dxy_chr6<-read.table("datafiles/chr6_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")
#Miss20Maf0.01
table_dxy_chr12<-read.table("datafiles/chr12_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")
#Miss20Maf0.01
table_dxy_chr17<-read.table("datafiles/chr17_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")
#Miss20Maf0.01
table_dxy_chr23<-read.table("datafiles/chr23_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")
#Miss20Maf0.01
table_dxy_ATL_PAC<-read.table("datafiles/atlantic_vs_pacific_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_dxy.txt",header=T,sep="\t")

table_dxy_chr6$mid<-(table_dxy_chr6$window_pos_1 + table_dxy_chr6$window_pos_2 )/2
table_dxy_chr12$mid<-(table_dxy_chr12$window_pos_1 + table_dxy_chr12$window_pos_2 )/2
table_dxy_chr17$mid<-(table_dxy_chr17$window_pos_1 + table_dxy_chr17$window_pos_2 )/2
table_dxy_chr23$mid<-(table_dxy_chr23$window_pos_1 + table_dxy_chr23$window_pos_2 )/2
table_dxy_ATL_PAC$mid<-(table_dxy_ATL_PAC$window_pos_1 + table_dxy_ATL_PAC$window_pos_2 )/2

table_dxy_chr6<-table_dxy_chr6 %>% filter(no_sites>=10000)
table_dxy_chr12<-table_dxy_chr12 %>% filter(no_sites>=10000)
table_dxy_chr17<-table_dxy_chr17 %>% filter(no_sites>=10000)
table_dxy_chr23<-table_dxy_chr23 %>% filter(no_sites>=10000)
table_dxy_ATL_PAC<-table_dxy_ATL_PAC %>% filter(no_sites>=10000)

#### Chr6 ####
s=22282765
e=24868581
chr="chr6"
POP1=c("CS_chr6")
POP2=c("BS_chr6")
XMinLim=s/1000000-2
XMaxLim=e/1000000+2

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

quant_0.99<-table_dxy_chr6 %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  summarise(quantile_CS_F=quantile(avg_dxy, 0.99, na.rm=T))

quant_0.95<-table_dxy_chr6 %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  summarise(quantile_CS_F=quantile(avg_dxy, 0.95, na.rm=T))


chr6_dxy<-table_dxy_chr6 %>%
  filter(chromosome==chr) %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.02),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(data=(table_dxy_ATL_PAC %>% filter(chromosome==chr)), aes(x=mid/1000000,y=avg_dxy), color="deepskyblue", size=1)+
  geom_line(aes(x=mid/1000000,y=avg_dxy),color="black",size=1)+
  
  geom_hline(yintercept = quant_0.95$quantile_CS_F,color="black",linetype="dashed")+
  xlim(XMinLim,XMaxLim)+
  ylim(0,0.02)+
  #ggtitle("pixy")+
  ylab("Dxy") +
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  my_thm_bottom_plot

#### Chr12 ####
s=17826318
e=25603093
chr="chr12"
POP1=c("CS_chr12")
POP2=c("BS_chr12")
XMinLim=s/1000000-2
XMaxLim=e/1000000+2

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

quant_0.95<-table_dxy_chr12 %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  summarise(quantile_CS_F=quantile(avg_dxy, 0.95, na.rm=T))

chr12_dxy<-table_dxy_chr12 %>%
  filter(chromosome==chr) %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.02),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  
  geom_line(data=(table_dxy_ATL_PAC %>% filter(chromosome==chr)), aes(x=mid/1000000,y=avg_dxy), color="deepskyblue", size=1)+
  geom_line(aes(x=mid/1000000,y=avg_dxy),color="black",size=1)+
  
  geom_hline(yintercept = quant_0.95$quantile_CS_F,color="black",linetype="dashed")+
  xlim(XMinLim,XMaxLim)+
  ylim(0,0.02)+
  #ggtitle("pixy")+
  ylab("Dxy") +
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  my_thm_bottom_plot

#### Chr17 ####
s=25805445
e=27568511
chr="chr17"
POP1=c("CS_chr17")
POP2=c("BS_chr17")
XMinLim=s/1000000-2
XMaxLim=e/1000000

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

quant_0.95<-table_dxy_chr17 %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  summarise(quantile_CS_F=quantile(avg_dxy, 0.95, na.rm=T))

chr17_dxy<-table_dxy_chr17 %>%
  filter(chromosome==chr) %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.02),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  # Plot Atlantic vs Pacific 
  geom_line(data=(table_dxy_ATL_PAC %>% filter(chromosome==chr)), aes(x=mid/1000000,y=avg_dxy), color="deepskyblue", size=1)+
  # Plot N vs S
  geom_line(aes(x=mid/1000000,y=avg_dxy),color="black",size=1)+
  
  geom_hline(yintercept = quant_0.95$quantile_CS_F,color="black",linetype="dashed")+
  xlim(XMinLim,XMaxLim)+
  ylim(0,0.02)+
  #ggtitle("pixy")+
  ylab("Dxy") +
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  my_thm_bottom_plot

#### Chr23 ####
s=16226443
e=17604273
chr="chr23"
POP1=c("CS_chr23")
POP2=c("F_chr23")
XMinLim=s/1000000-2
XMaxLim=e/1000000+2


candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

quant_0.95<-table_dxy_chr23 %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  summarise(quantile_CS_F=quantile(avg_dxy, 0.95, na.rm=T))


chr23_dxy<-table_dxy_chr23 %>%
  filter(chromosome==chr) %>%
  filter(pop1==POP1 & pop2==POP2) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.02),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(data=(table_dxy_ATL_PAC %>% filter(chromosome==chr)), aes(x=mid/1000000,y=avg_dxy), color="deepskyblue", size=1)+
  geom_line(aes(x=mid/1000000,y=avg_dxy),color="black",size=1)+
  
  geom_hline(yintercept = quant_0.95$quantile_CS_F,color="black",linetype="dashed")+
  xlim(XMinLim,XMaxLim)+
  ylim(0,0.02)+
  #ggtitle("pixy")+
  ylab("Dxy") +
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  my_thm_bottom_plot

all_plots_dxy<-grid.arrange(chr6_dxy, chr12_dxy, chr17_dxy, chr23_dxy)

ggsave(all_plots_dxy, file="figures/all_plots_dxy_comparison_Atlantic_Pacific_herring.pdf", unit="mm", width=160, height =100)

## Supplementary Figure 10B. Nucleotide diversity of all Atlantic herring populations ####

# Note: In the end, data from these chromosome specific files were not used in the plot, but since I didn't modify the code 
# the files are still necessary for the plots bellow. 

# Miss 20 Maf 0.01
table_pi_chr6<-read.table("datafiles/chr6_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)
# Miss 20 Maf 0.01
table_pi_chr12<-read.table("datafiles/chr12_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)
# Miss 20 Maf 0.01
table_pi_chr17<-read.table("datafiles/chr17_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)
# Miss 20 Maf 0.01
table_pi_chr23<-read.table("datafiles/chr23_inversion_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt", header=T)
# Miss 20 Maf 0.01
table_pi_ATL_PAC<-read.table("datafiles/atlantic_vs_pacific_popgen_stats.repetition_miss20maf0.01.wg.20kb.popgenpixy.out_pi.txt",header=T,sep="\t")

table_pi_chr6$mid<-(table_pi_chr6$window_pos_1 + table_pi_chr6$window_pos_2 )/2
table_pi_chr12$mid<-(table_pi_chr12$window_pos_1 + table_pi_chr12$window_pos_2 )/2
table_pi_chr17$mid<-(table_pi_chr17$window_pos_1 + table_pi_chr17$window_pos_2 )/2
table_pi_chr23$mid<-(table_pi_chr23$window_pos_1 + table_pi_chr23$window_pos_2 )/2
table_pi_ATL_PAC$mid<-(table_pi_ATL_PAC$window_pos_1 + table_pi_ATL_PAC$window_pos_2 )/2

table_pi_chr6 <-table_pi_chr6 %>% filter(no_sites>=10000)
table_pi_chr12 <-table_pi_chr12 %>% filter(no_sites>=10000)
table_pi_chr17 <-table_pi_chr17 %>% filter(no_sites>=10000)
table_pi_chr23 <-table_pi_chr23 %>% filter(no_sites>=10000)
table_pi_ATL_PAC <-table_pi_ATL_PAC %>% filter(no_sites>=10000)


#### Chr6 ####
s=22282765
e=24868581
chr="chr6"
POP1=c("CS_chr6")
POP2=c("BS_chr6")
XMinLim=s/1000000-2
XMaxLim=e/1000000+2

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

quant_0.99_pi<-table_pi_chr6 %>%
  filter(pop!="Vancouver") %>%
  summarise(quantile_CS_F=quantile(avg_pi, 0.99, na.rm=T))

pi_chr6<-table_pi_chr6 %>%
  filter(chromosome==chr) %>%
  filter(pop!="Vancouver") %>%
  filter(pop==POP1 | pop==POP2) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.01),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(data=(table_pi_ATL_PAC %>% filter(pop=="Atlantic" & chromosome==chr)), aes(x=mid/1000000, y=avg_pi), color="gray25", size=1)+
  #geom_line(aes(x=mid/1000000, y=avg_pi, color=pop), size=1, alpha=0.5)+
  #scale_color_manual(values=c("#1F78B4", "#D95F02"))+
  #geom_hline(yintercept = quant_0.99_pi$quantile_CS_F, color="black",linetype="dashed")+
  #xlim(XMinLim,XMaxLim)+
  ylim(0,0.01)+
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  ylab("Pi")+
  my_thm_bottom_plot

#### Chr12 ####
s=17826318
e=25603093
chr="chr12"
POP1=c("CS_chr12")
POP2=c("BS_chr12")
XMinLim=s/1000000-2
XMaxLim=e/1000000+2

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

quant_0.99_pi<-table_pi_chr12 %>%
  filter(pop!="Vancouver") %>%
  summarise(quantile_CS_F=quantile(avg_pi, 0.99, na.rm=T))

pi_chr12<-table_pi_chr12 %>%
  filter(chromosome==chr) %>%
  filter(pop!="Vancouver") %>%
  filter(pop==POP1 | pop==POP2) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.01),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(data=(table_pi_ATL_PAC %>% filter(pop=="Atlantic" & chromosome==chr)), aes(x=mid/1000000, y=avg_pi), color="gray25", size=1)+
  #geom_line(aes(x=mid/1000000, y=avg_pi, color=pop), size=1, alpha=0.5)+
  #scale_color_manual(values=c("#1F78B4", "#D95F02"))+
  #geom_hline(yintercept = quant_0.99_pi$quantile_CS_F, color="black",linetype="dashed")+
  #xlim(XMinLim,XMaxLim)+
  ylim(0,0.01)+
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  ylab("Pi")+
  my_thm_bottom_plot


#### Chr17 ####
s=25805445
e=27568511
chr="chr17"
POP1=c("CS_chr17")
POP2=c("BS_chr17")
XMinLim=s/1000000-2
XMaxLim=e/1000000

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

quant_0.99_pi<-table_pi_chr17 %>%
  filter(pop!="Vancouver") %>%
  summarise(quantile_CS_F=quantile(avg_pi, 0.99, na.rm=T))

pi_chr17<-table_pi_chr17 %>%
  filter(chromosome==chr) %>%
  filter(pop!="Vancouver") %>%
  filter(pop==POP1 | pop==POP2) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.01),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(data=(table_pi_ATL_PAC %>% filter(pop=="Atlantic" & chromosome==chr)), aes(x=mid/1000000, y=avg_pi), color="gray25", size=1)+
  #geom_line(aes(x=mid/1000000, y=avg_pi, color=pop), size=1, alpha=0.5)+
  #scale_color_manual(values=c("#1F78B4", "#D95F02"))+
  #geom_hline(yintercept = quant_0.99_pi$quantile_CS_F, color="black",linetype="dashed")+
  #xlim(XMinLim,XMaxLim)+
  ylim(0,0.01)+
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  ylab("Pi")+
  my_thm_bottom_plot

#### Chr23 ####
s=16226443
e=17604273
chr="chr23"
POP1=c("CS_chr23")
POP2=c("F_chr23")
XMinLim=s/1000000-2
XMaxLim=e/1000000+2

candidates<-data.frame(matrix(c(s/1000000,e/1000000),ncol=2,nrow=1))

quant_0.99_pi<-table_pi_chr23 %>%
  filter(pop!="Vancouver") %>%
  summarise(quantile_CS_F=quantile(avg_pi, 0.99, na.rm=T))

pi_chr23<-table_pi_chr23 %>%
  filter(chromosome==chr) %>%
  filter(pop!="Vancouver") %>%
  filter(pop==POP1 | pop==POP2) %>%
  ggplot()+
  geom_rect(data=candidates,mapping=aes(xmin=X1,xmax=X2,ymin=0,ymax=0.01),
            fill="gray",color="NA",alpha=0.5,inherit.aes = FALSE)+
  geom_line(data=(table_pi_ATL_PAC %>% filter(pop=="Atlantic" & chromosome==chr)), aes(x=mid/1000000, y=avg_pi), color="gray25", size=1)+
  #geom_line(aes(x=mid/1000000, y=avg_pi, color=pop), size=1, alpha=0.5)+
  #scale_color_manual(values=c("#1F78B4", "#D95F02"))+
  #geom_hline(yintercept = quant_0.99_pi$quantile_CS_F, color="black",linetype="dashed")+
  #xlim(XMinLim,XMaxLim)+
  ylim(0,0.01)+
  xlab(paste0("Chromosome ", str_remove(chr, "chr"), " (Mb)")) +
  ylab("Pi")+
  my_thm_bottom_plot

all_pi_plots<-grid.arrange(pi_chr6, pi_chr12, pi_chr17, pi_chr23, nrow=4)

ggsave(all_pi_plots, file="figures/chromosome_wide_pi_plots_comparison_with_Atlantic.pdf", units="mm", width=160, height = 100)
ggsave(all_pi_plots, file="figures/chromosome_wide_pi_plots_Atlantic.pdf", units="mm", width=160, height = 100)


# Figure 5E - Genome-wide and Northern and Southern haplotype nucleotide diversity ####

# Boxplots
table_pi_inversion_chr6$type<-"inversion_chr6"
table_pi_inversion_chr12$type<-"inversion_chr12"
table_pi_inversion_chr17$type<-"inversion_chr17"
table_pi_inversion_chr23$type<-"inversion_chr23"

table_pi_inversion_chr6$ancestral<-NA
table_pi_inversion_chr12$ancestral<-NA
table_pi_inversion_chr17$ancestral<-NA
table_pi_inversion_chr23$ancestral<-NA

table_pi_inversion_chr6[table_pi_inversion_chr6$pop=="CS_chr6",]$ancestral<-"ancestral"
table_pi_inversion_chr6[table_pi_inversion_chr6$pop=="BS_chr6",]$ancestral<-"derived"

table_pi_inversion_chr12[table_pi_inversion_chr12$pop=="CS_chr12",]$ancestral<-"ancestral"
table_pi_inversion_chr12[table_pi_inversion_chr12$pop=="BS_chr12",]$ancestral<-"derived"

table_pi_inversion_chr17[table_pi_inversion_chr17$pop=="CS_chr17",]$ancestral<-"derived"
table_pi_inversion_chr17[table_pi_inversion_chr17$pop=="BS_chr17",]$ancestral<-"ancestral"

table_pi_inversion_chr23[table_pi_inversion_chr23$pop=="CS_chr23",]$ancestral<-"derived"
table_pi_inversion_chr23[table_pi_inversion_chr23$pop=="F_chr23",]$ancestral<-"ancestral"

table_pi_inversion_chr6$haplotype<-NA
table_pi_inversion_chr12$haplotype<-NA
table_pi_inversion_chr17$haplotype<-NA
table_pi_inversion_chr23$haplotype<-NA

table_pi_inversion_chr6[table_pi_inversion_chr6$pop=="CS_chr6",]$haplotype<-"South"
table_pi_inversion_chr6[table_pi_inversion_chr6$pop=="BS_chr6",]$haplotype<-"North"

table_pi_inversion_chr12[table_pi_inversion_chr12$pop=="CS_chr12",]$haplotype<-"South"
table_pi_inversion_chr12[table_pi_inversion_chr12$pop=="BS_chr12",]$haplotype<-"North"

table_pi_inversion_chr17[table_pi_inversion_chr17$pop=="CS_chr17",]$haplotype<-"South"
table_pi_inversion_chr17[table_pi_inversion_chr17$pop=="BS_chr17",]$haplotype<-"North"

table_pi_inversion_chr23[table_pi_inversion_chr23$pop=="CS_chr23",]$haplotype<-"South"
table_pi_inversion_chr23[table_pi_inversion_chr23$pop=="F_chr23",]$haplotype<-"North"


inversion<-rbind(table_pi_inversion_chr6, table_pi_inversion_chr12, table_pi_inversion_chr17, table_pi_inversion_chr23)

table_pi_ATL$type<-"genome_wide"
table_pi_ATL$ancestral<-"genome_wide"
table_pi_ATL$haplotype<-"genome_wide"

to_plot_pi_boxplot<-rbind(table_pi_ATL, inversion)

to_plot_pi_boxplot$pop<-factor(to_plot_pi_boxplot$pop, 
                               levels=c("Atlantic", 
                                        "BS_chr6", "CS_chr6", 
                                        "BS_chr12", "CS_chr12", 
                                        "BS_chr17", "CS_chr17",
                                        "F_chr23", "CS_chr23"))

to_plot_stats$xx<-c("CS_chr6", "BS_chr6", "CS_chr12", "BS_chr12", "CS_chr17", "BS_chr17", "CS_chr23", "F_chr23", "Atlantic")
to_plot_stats$haplotype<-c("South", "BS_chr6", "CS_chr12", "BS_chr12", "CS_chr17", "BS_chr17", "CS_chr23", "F_chr23", "Atlantic")

geom_boxplot_pi<-ggplot(to_plot_pi_boxplot, aes(x=pop, y=avg_pi*100)) +
  geom_boxplot(aes(fill=haplotype), lwd=0.2)+
  geom_point(data=to_plot_stats, aes(y=1.5, x=xx, shape=ancestral, color=ancestral), fill="black", size=2)+
  geom_text(data=to_plot_stats, aes(xx, values, label =sprintf("%0.3f", round(values, digits = 3))), vjust = -0.5)+
  scale_fill_manual(values=c("white", "#1F78B4", "#D95F02"))+
  scale_color_manual(values=c("darkgreen", "#7030A0", "black"))+
  theme_classic()+
  ylab("Nucleotide diversity (Pi %)")+
  my_thm3

ggsave(geom_boxplot_pi, filename="figures/boxplot_distribution_pi_miss20_maf0.01_genome_vs_inversions.pdf", height = 1.5, width=8, dpi=300) 

#### Wilcoxon TESTS ####

t.test(avg_pi~pop, data=table_pi_inversion_chr6)
#W = 6907, p-value = 0.07488
wilcox.test(avg_pi~pop, data=table_pi_inversion_chr12)
#W = 49296, p-value = 0.0003782
wilcox.test(avg_pi~pop, data=table_pi_inversion_chr17)
#W = 960, p-value = 5.38e-10
wilcox.test(avg_pi~pop, data=table_pi_inversion_chr23)
#W = 1871, p-value = 0.003235

wilcox.test(table_pi_inversion_chr6[table_pi_inversion_chr6$pop=="BS_chr6",]$avg_pi, table_pi_ATL$avg_pi)
#W = 1398579, p-value = 5.587e-06
wilcox.test(table_pi_inversion_chr6[table_pi_inversion_chr6$pop=="CS_chr6",]$avg_pi, table_pi_ATL$avg_pi)
#W = 1520659, p-value = 0.001177
wilcox.test(table_pi_inversion_chr12[table_pi_inversion_chr12$pop=="BS_chr12",]$avg_pi, table_pi_ATL$avg_pi)
#W = 5859412, p-value = 6.58e-09
wilcox.test(table_pi_inversion_chr12[table_pi_inversion_chr12$pop=="CS_chr12",]$avg_pi, table_pi_ATL$avg_pi)
#W = 6195374, p-value = 1.724e-15
wilcox.test(table_pi_inversion_chr17[table_pi_inversion_chr17$pop=="BS_chr17",]$avg_pi, table_pi_ATL$avg_pi)
#W = 1614596, p-value < 2.2e-16
wilcox.test(table_pi_inversion_chr17[table_pi_inversion_chr17$pop=="CS_chr17",]$avg_pi, table_pi_ATL$avg_pi)
#W = 1818067, p-value < 2.2e-16
wilcox.test(table_pi_inversion_chr23[table_pi_inversion_chr23$pop=="F_chr23",]$avg_pi, table_pi_ATL$avg_pi)
#W = 905565, p-value = 0.02416
wilcox.test(table_pi_inversion_chr23[table_pi_inversion_chr23$pop=="CS_chr23",]$avg_pi, table_pi_ATL$avg_pi)
#W = 1074759, p-value = 4.917e-07

t.test(table_pi_inversion_chr6[table_pi_inversion_chr6$pop=="BS_chr6",]$avg_pi, table_pi_ATL$avg_pi)
#t = -16.815, df = 135.69, p-value < 2.2e-16
t.test(table_pi_inversion_chr6[table_pi_inversion_chr6$pop=="CS_chr6",]$avg_pi, table_pi_ATL$avg_pi)
#t = -14.508, df = 135.39, p-value < 2.2e-16
t.test(table_pi_inversion_chr12[table_pi_inversion_chr12$pop=="BS_chr12",]$avg_pi, table_pi_ATL$avg_pi)
#t = -0.80307, df = 376.26, p-value = 0.4224
t.test(table_pi_inversion_chr12[table_pi_inversion_chr12$pop=="CS_chr12",]$avg_pi, table_pi_ATL$avg_pi)
#t = 4.8032, df = 357.79, p-value = 2.301e-06
t.test(table_pi_inversion_chr17[table_pi_inversion_chr17$pop=="BS_chr17",]$avg_pi, table_pi_ATL$avg_pi)
#t = 10.096, df = 69.754, p-value = 2.79e-15
t.test(table_pi_inversion_chr17[table_pi_inversion_chr17$pop=="CS_chr17",]$avg_pi, table_pi_ATL$avg_pi)
#t = 18.65, df = 69.714, p-value < 2.2e-16
t.test(table_pi_inversion_chr23[table_pi_inversion_chr23$pop=="F_chr23",]$avg_pi, table_pi_ATL$avg_pi)
#t = 0.23729, df = 52.584, p-value = 0.8134
t.test(table_pi_inversion_chr23[table_pi_inversion_chr23$pop=="CS_chr23",]$avg_pi, table_pi_ATL$avg_pi)
#t = 4.2399, df = 52.275, p-value = 9.135e-05


# Supplementary Figure 15 - Genome wide Fst and Nm ####

table_fst<-read.table("datafiles/all_herring_populations.wg.20kb.popgenpixy.out_fst.txt", header=T)

head(table_fst)

# Exclude inversions regions

chr6<-table_fst[table_fst$chromosome=="chr6" & table_fst$window_pos_1>=22282765 & table_fst$window_pos_2 <=24868581,]
chr12<-table_fst[table_fst$chromosome=="chr12" & table_fst$window_pos_1>=17826318 & table_fst$window_pos_2 <=25603093,]
chr17<-table_fst[table_fst$chromosome=="chr17" & table_fst$window_pos_1>=25805445 & table_fst$window_pos_2 <=27568511,]
chr23<-table_fst[table_fst$chromosome=="chr23" & table_fst$window_pos_1>=16226443 & table_fst$window_pos_2 <=17604273,]

inversions<-rbind(chr6, chr12, chr17, chr23)
to_exclude<-rownames(inversions)
all_rows<-rownames(table_fst)
all_rows %in% to_exclude

table_fst_without_inversions<-table_fst[!(all_rows %in% to_exclude),]

fst_all_pops<-table_fst_without_inversions %>% 
  group_by(pop1, pop2) %>%
  filter(avg_wc_fst >= 0) %>%
  summarise(mean(avg_wc_fst, na.rm=T))

# CS vs BS 
# 0.02221297
# Using Hoekstra 2006

#Nm=[(1/FST)−1]/4

((1/0.02221297)-1)/4
11.00469

# for all pops, calculate Nm:
fst_all_pops$pop1<-factor(fst_all_pops$pop1, levels=c("Baltic_Spring", "Baltic_Autumn", "NorthSea", "UK", "Norway_Spring", "Canada_Spring", "Canada_Autumn"))
fst_all_pops$pop2<-factor(fst_all_pops$pop2, levels=c("Baltic_Spring", "Baltic_Autumn", "NorthSea", "UK", "Norway_Spring", "Canada_Spring", "Canada_Autumn"))

fst_all_pops$Nm <- ((1/fst_all_pops$`mean(avg_wc_fst, na.rm = T)`)-1)/4

fst<-ggplot(fst_all_pops, aes(x=pop1, pop2, fill=`mean(avg_wc_fst, na.rm = T)`))+
  geom_tile()+
  geom_text(aes(label = round(`mean(avg_wc_fst, na.rm = T)`, 3)), color = "white", size = 3)+
  scale_fill_gradient(low = "blue", high = "red", name="Fst")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90, size=10),
        axis.text.y=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

Nm<-ggplot(fst_all_pops, aes(x=pop1, pop2, fill=Nm))+
  geom_tile()+
  geom_text(aes(label = round(Nm, 2)), color = "white", size = 3)+
  scale_fill_gradient(low = "blue", high = "red", name="Nm")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90, size=10),
        axis.text.y=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

Nm_fst<-gridExtra::grid.arrange(fst, Nm, nrow=1)

ggsave(plot=Nm_fst, filename="figures/pairwise_Nm_Fst_AtlHerring_populations.pdf", width=10, height = 5)

