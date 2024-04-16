# Run circular visiualization in R

install.packages("circlize")
library(circlize)
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(reshape2)
library(gridExtra)
library(ggrastr)

session<-sessionInfo()

save.image("~/Documents/Postdoc/Repositories/Atlherr_inversion_evol_RData/synteny.RData")

satsuma<- read.table("alignments/satsuma_summary.chained.out_sorted", stringsAsFactors=F, sep = "\t", comment.char = "")
colnames(satsuma)<-c("q_scaffold", "q_start", "q_end", "t_scaffold", "t_start", "t_end", "NN", "strand")

sprat_index<-read.table("data/GCA_963457725.1_fSprSpr1.1_genomic.fna.fai")
herring_index<-read.table("data/Ch_v2.0.2.fasta.fai")

# select the scaffolds with the most mappings to herring:
tt<-satsuma %>% group_by(q_scaffold) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  filter(n>1e5)

print(tt, n=200)

# subset the genomes to the largest scaffolds/ chromosomes:
herring_index_chr<-head(herring_index, n=26)

tt$q_scaffold_clean<- str_split_fixed(tt$q_scaffold, "_", 2)[,1]
tt$q_chromosome<- str_split_fixed(tt$q_scaffold, "_", 7)[,7]
sprat_index_chr<- sprat_index %>% filter(V1 %in% tt$q_scaffold_clean)
sprat_index_chr$chr<-tt$q_chromosome

herring_df<-data.frame(name=herring_index_chr$V1, start=1, end=herring_index_chr$V2)
sprat_df<-data.frame(scaffold=sprat_index_chr$V1, start=1, end=sprat_index_chr$V2)
sprat_df$name<-sprat_index_chr$chr
sprat_df<-sprat_df[,c(4,2,3,1)]

# combine genomes:

combined_genomes <- rbind(herring_df, sprat_df[,c(1,2,3)])

chromosome.index = c(paste0("chr", c(1:26)), 
                     rev(paste0("scaf", c(1:20))))

satsuma$sprat_scaf<- str_split_fixed(satsuma$q_scaffold, "_", 2)[,1]

# filter regions that are going to be in the plot:
satsuma_filter<-satsuma %>% 
  filter(t_scaffold %in% herring_df$name) %>%
  filter(sprat_scaf %in% sprat_df$scaffold)

# calculate a proportion of each herring chromosome covered by each sprat scaffold

satsuma_filter_chrLen<-merge(satsuma_filter, herring_index, by.x="t_scaffold", by.y="V1")

satsuma_filter_chrLen$t_proportion<- (satsuma_filter_chrLen$t_end - satsuma_filter_chrLen$t_start)/ satsuma_filter_chrLen$V2

df_prop<-satsuma_filter_chrLen %>% 
  group_by(t_scaffold, sprat_scaf) %>%
  summarise(value=sum(t_proportion)) %>%
  arrange(desc(value))

print(satsuma_filter_chrLen %>% 
        group_by(t_scaffold) %>%
        summarise(value=sum(t_proportion)), n=26)

colnames(df_prop)<-c("from", "to", "value") 

df_prop<-merge(df_prop, sprat_df, by.x="to", by.y="scaffold")
df_prop<-data.frame(from=df_prop$from, to=df_prop$name, value=df_prop$value)

df_prop$value<-df_prop$value*100
df_prop$to<-paste0("Cpr",df_prop$to)
df_prop$from<-str_replace(df_prop$from, "c", toupper)

chromosome.index = c(paste0("Cpr", c(20:1)),
                     rev(paste0("Chr", c(26:1))))

chordDiagram(df_prop, order=chromosome.index, annotationTrack = c("grid", "axis"))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1]+2, CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.5)
}, bg.border = NA) # here set bg.border to NA is important

circos.clear()

# Cool, but now for the three chromosomes of interest...

satsuma_filter_chr<-satsuma %>% filter(t_scaffold %in% c("chr6", "chr12", "chr17", "chr23"))

# what if we do a dotplot:
satsuma_filter_chr$t_mid<-(satsuma_filter_chr$t_end+satsuma_filter_chr$t_start)/2
satsuma_filter_chr$q_mid<-(satsuma_filter_chr$q_end+satsuma_filter_chr$q_start)/2

chr6<-satsuma_filter_chr %>% 
  filter(t_scaffold=="chr6") %>%
  ggplot()+
  geom_rect(mapping=aes(xmin=22282765/1e6,xmax=24868581/1e6,ymin=0, ymax=max(q_end)/1e6), color=NA, fill="gray", alpha=0.5)+
  rasterise(geom_point(aes(x=t_mid/1e6, y=q_mid/1e6, color=strand)), dpi=300)+
  theme_classic()+
  xlab("Position Chromosome 6 Herring (Mb)")+
  ylab("Position Sprat (Mb)")+
  theme(legend.position="none")

chr12<-satsuma_filter_chr %>% 
  filter(t_scaffold=="chr12") %>%
  
  ggplot()+
  geom_rect(mapping=aes(xmin=17826318/1e6,xmax=25603093/1e6,ymin=0, ymax=max(q_end)/1e6), color=NA, fill="gray", alpha=0.5)+
  rasterise(geom_point(aes(x=t_mid/1e6, y=q_mid/1e6, color=strand)), dpi=300)+
  theme_classic()+
  xlab("Position Chromosome 12 Herring (Mb)")+
  ylab("Position Sprat (Mb)")+
  theme(legend.position="none")

chr17<-satsuma_filter_chr %>% 
  filter(t_scaffold=="chr17") %>%
  ggplot()+
  rasterise(geom_rect(mapping=aes(xmin=25805445/1e6,xmax=27568511/1e6,ymin=0, ymax=max(q_end)/1e6), color=NA, fill="gray", alpha=0.5), dpi=300)+
  geom_point(aes(x=t_mid, y=q_mid, color=strand))+
  theme_classic()+
  xlab("Position Chromosome 17 Herring (Mb)")+
  ylab("Position Sprat (MB)")+
  theme(legend.position="none")

chr23<-satsuma_filter_chr %>% 
  filter(t_scaffold=="chr23") %>%
  ggplot()+
  geom_rect(mapping=aes(xmin=16226443/1e6,xmax=17604273/1e6,ymin=0, ymax=max(q_end)/1e6), color=NA, fill="gray", alpha=0.5)+
  rasterise(geom_point(aes(x=t_mid/1e6, y=q_mid/1e6, color=strand)), dpi=300)+
  theme_classic()+
  xlab("Position Chromosome 23 Herring (Mb)")+
  ylab("Position Sprat (Mb)")+
  theme(legend.position="none")

inversion_plots<-grid.arrange(chr6, chr12, chr17, chr23, nrow=2)

ggsave(inversion_plots, filename="figures/dot_plot_chr6_chr12_chr17_chr23.pdf", unit="mm", width=180, height=180)
ggsave(inversion_plots, filename="figures/dot_plot_chr6_chr12_chr17_chr23.png", unit="mm", width=180, height=180, dpi = 300)


