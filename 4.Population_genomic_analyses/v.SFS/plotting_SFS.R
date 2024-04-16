# Code to Plot SFS from the outputs of vcftools

# Load libraries ####
require(tidyverse)

# Set working directory ####
setwd("~/Documents/Postdoc/Project_Herring/Inversion_Project/Mutation_load/SFS/codon_positions_analysis_v02_2023-03-14")

# Read input files ####
nonsyn_files<-list.files(path = "data/", pattern="_nonsyn.mod.frq")
syn_files<-list.files(path = "data/", pattern="_syn.mod.frq")

output_names<-str_remove(nonsyn_files, pattern = "_nonsyn.mod.frq")
chr<-str_split_fixed(nonsyn_files, pattern="_", n=4)[,3]

# This loop will go through the files and plot the SFS for each inversion and haplotype

for(i in 1:length(nonsyn_files)){
  # Let's try this as counts:
  nonsyn_input<-read.table(paste0("data/",nonsyn_files[i]), skip=1)
  syn_input<-read.table(paste0("data/",syn_files[i]), skip=1)
  
  # Change the column headers:
  colnames(nonsyn_input)<-c("chr", "pos", "no_alleles", "no_chr", "Anc", "freqAnc", "Dev", "freqDev")
  colnames(syn_input)<-c("chr", "pos", "no_alleles", "no_chr", "Anc", "freqAnc", "Dev", "freqDev")
  
  # FIlter SFS to remove fixed sites. 
  
  nonsyn_input_polym<-nonsyn_input %>% 
    filter(freqDev!=1) %>% 
    filter(freqAnc!=1)
  
  syn_input_polym<-syn_input %>% 
    filter(freqDev!=1) %>% 
    filter(freqAnc!=1)
  
  # Let's try a tidyverse approach
  
  nonsyn_input_polym_bins<-nonsyn_input_polym %>%
    mutate(bin = cut(freqDev, seq(0, 1, 0.1), right = FALSE, include.lowest=TRUE))
  
  
  syn_input_polym_bins<-syn_input_polym %>%
    mutate(bin = cut(freqDev, seq(0, 1, 0.1), right = FALSE, include.lowest=TRUE))
  
  # Calculate counts per bin, and add a column with count type
  nonsyn_input_polym_counts<-nonsyn_input_polym_bins %>% group_by(bin) %>% summarise(n=n()) %>% add_column(type="NonSyn")
  
  syn_input_polym_counts<-syn_input_polym_bins %>% group_by(bin) %>% summarise(n=n()) %>%        add_column(type="Syn")
  
  to_plot_counts <- rbind(nonsyn_input_polym_counts, syn_input_polym_counts)
  
  # Plot
  SFS<-ggplot(data = to_plot_counts, mapping = aes(y=n, x=bin, fill=type)) + 
    geom_bar(color="white",alpha=0.7, stat = "identity", position="dodge") +
    labs(x='derived allele frequency', y="SNP counts") +
    ggtitle(paste0(output_names[i]," \nNo. NONSYN SNPs: ", nrow(nonsyn_input_polym), " \nNo. SYN SNPs: ", nrow(syn_input_polym)))+
    theme_minimal() 
  
  ggsave(plot = SFS, filename=paste0("figures/",output_names[i],".coding.all.chromosome.SFS.pdf"), height = 4, width=12)
  ggsave(plot = SFS, filename=paste0("figures/",output_names[i],".coding.all.chromosome.SFS.pdf"), height = 4, width=12)
}

# Plotting a joint plot for Figure 6 ####
chr6_inv<-"chr6"
start_chr6<-22282765	
end_chr6<-24868581

chr12_inv<-"chr12"
start_chr12<-17826318
end_chr12<-25603093

chr17_inv<-"chr17"	
start_chr17<-25805445	
end_chr17<-27568511

chr23_inv<-"chr23"
start_chr23<-16226443
end_chr23<-17604273

inversion_coordinates<-data.frame(CHR=c(6, 12, 17, 23), 
           START=c(start_chr6, start_chr12, start_chr17, start_chr23),
           END=c(end_chr6, end_chr12, end_chr17, end_chr23))

plot_list = list()
for(i in 1:length(nonsyn_files)){
  #
  # Let's try this as counts:
  nonsyn_input<-read.table(paste0("data/",nonsyn_files[i]), skip=1)
  syn_input<-read.table(paste0("data/",syn_files[i]), skip=1)
  
  # Change the column headers:
  colnames(nonsyn_input)<-c("chr", "pos", "no_alleles", "no_chr", "Anc", "freqAnc", "Dev", "freqDev")
  colnames(syn_input)<-c("chr", "pos", "no_alleles", "no_chr", "Anc", "freqAnc", "Dev", "freqDev")
  
  # Filter input table to just look at inversion coordinates
  inversion<-str_remove(str_split_fixed(nonsyn_files[i], "_", 4)[,3], "chr")
  df_coord<-inversion_coordinates %>% filter(CHR==inversion)
  
  nonsyn_input<-nonsyn_input %>% filter(chr==df_coord[1,"CHR"] & pos>=df_coord[1,"START"] & pos<=df_coord[1,"END"])
  syn_input<-syn_input %>% filter(chr==df_coord[1,"CHR"] & pos>=df_coord[1,"START"] & pos<=df_coord[1,"END"])
  
  # FIlter SFS to remove fixed sites. 
  
  nonsyn_input_polym<-nonsyn_input %>% 
    filter(freqDev!=1) %>% 
    filter(freqAnc!=1)
  
  syn_input_polym<-syn_input %>% 
    filter(freqDev!=1) %>% 
    filter(freqAnc!=1)
  
  # Let's try a tidyverse approach
  
  nonsyn_input_polym_bins<-nonsyn_input_polym %>%
    mutate(bin = cut(freqDev, seq(0, 1, 0.1), right = FALSE, include.lowest=TRUE))
  
  
  syn_input_polym_bins<-syn_input_polym %>%
    mutate(bin = cut(freqDev, seq(0, 1, 0.1), right = FALSE, include.lowest=TRUE))
  
  # Calculate counts per bin, and add a column with count type
  nonsyn_input_polym_counts<-nonsyn_input_polym_bins %>% group_by(bin) %>% summarise(n=n()) %>% add_column(type="NonSyn")
  
  syn_input_polym_counts<-syn_input_polym_bins %>% group_by(bin) %>% summarise(n=n()) %>%        add_column(type="Syn")
  
  to_plot_counts <- rbind(nonsyn_input_polym_counts, syn_input_polym_counts)
  
  to_plot_counts$proportion<-to_plot_counts$n/sum(to_plot_counts$n)
  
  # Plot
  SFS<-ggplot(data = to_plot_counts, mapping = aes(y=proportion, x=bin, fill=type)) + 
    geom_bar(color="white",alpha=0.7, stat = "identity", position="dodge") +
    labs(x='derived allele frequency', y="SNP counts") +
    #ggtitle(paste0(output_names[i]," \nNo. NONSYN SNPs: ", nrow(nonsyn_input_polym), " \nNo. SYN SNPs: ", nrow(syn_input_polym)))+
    theme_classic()+
    ylim(0,0.5)+
    theme(title = element_text(size=7),
          legend.position = "none",
          axis.text.x = element_text(size=7, angle=90),
          axis.text.y = element_text(size=7),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  plot_list[[i]] = SFS
  
}

plot_all<-gridExtra::grid.arrange(plot_list[[4]], plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[8]], plot_list[[5]], plot_list[[6]], plot_list[[7]], nrow=2, ncol=4)

ggsave(plot_all, file="figures/all_SFS_plots.pdf", height =100, width=198, dpi=300, units = "mm")
ggsave(plot_all, file="figures/all_SFS_plots_with_annotations.pdf", height =100, width=198, dpi=300, units = "mm")

ggsave(plot_all, file="figures/all_SFS_plots_with_annotations_proportion_SNPs.pdf", height =100, width=198, dpi=300, units = "mm")


#### Plotting genome-wide SFS ####
setwd("~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Inversion_Project/Mutation_load/SFS/codon_positions_analysis_v02_2023-03-14/data/genome-wide")

nonsyn_files<-list.files(pattern="_nonsyn.mod.frq")
syn_files<-list.files(pattern="_syn.mod.frq")

#output_names<-str_remove(nonsyn_files, pattern = "_nonsyn.mod.frq")
#chr<-str_split_fixed(nonsyn_files, pattern="_", n=4)[,3]

# read all the tables:
nonsyn_tables<-lapply(nonsyn_files, read.table, skip=1)
syn_tables<-lapply(syn_files, read.table, skip=1)

# bind them by row:
nonsyn_wg<-bind_rows(nonsyn_tables)
syn_wg<-bind_rows(syn_tables)

# Change the column headers:
colnames(nonsyn_wg)<-c("chr", "pos", "no_alleles", "no_chr", "Anc", "freqAnc", "Dev", "freqDev")
colnames(syn_wg)<-c("chr", "pos", "no_alleles", "no_chr", "Anc", "freqAnc", "Dev", "freqDev")

# FIlter SFS to remove fixed sites. 

nonsyn_input_polym<-nonsyn_wg %>% 
  filter(freqDev!=1) %>% 
  filter(freqAnc!=1)

syn_input_polym<-syn_wg %>% 
  filter(freqDev!=1) %>% 
  filter(freqAnc!=1)

# Let's try a tidyverse approach

nonsyn_input_polym_bins<-nonsyn_input_polym %>%
  mutate(bin = cut(freqDev, seq(0, 1, 0.1), right = FALSE, include.lowest=TRUE))


syn_input_polym_bins<-syn_input_polym %>%
  mutate(bin = cut(freqDev, seq(0, 1, 0.1), right = FALSE, include.lowest=TRUE))

# Calculate counts per bin, and add a column with count type
nonsyn_input_polym_counts<-nonsyn_input_polym_bins %>% group_by(bin) %>% summarise(n=n()) %>% add_column(type="NonSyn")

syn_input_polym_counts<-syn_input_polym_bins %>% group_by(bin) %>% summarise(n=n()) %>%        add_column(type="Syn")

to_plot_counts <- rbind(nonsyn_input_polym_counts, syn_input_polym_counts)

to_plot_counts$proportion<-to_plot_counts$n/sum(to_plot_counts$n)

SFS<-ggplot(data = to_plot_counts, mapping = aes(y=proportion, x=bin, fill=type)) + 
  geom_bar(color="white",alpha=0.7, stat = "identity", position="dodge") +
  labs(x='derived allele frequency', y="SNP proportion") +
  #ggtitle(paste0(output_names[i]," \nNo. NONSYN SNPs: ", nrow(nonsyn_input_polym), " \nNo. SYN SNPs: ", nrow(syn_input_polym)))+
  theme_classic()+
  ylim(0,0.5)+
  theme(title = element_text(size=7),
        legend.position = "none",
        axis.text.x = element_text(size=7, angle=90),
        axis.text.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(SFS, file="../../figures/genome_wide_SFS_SNP_proportion.pdf", height =100, width=100, dpi=300, units = "mm")

# Session Info ####
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
#  [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.3     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1   
#  [9] ggplot2_3.4.4   tidyverse_2.0.0

# loaded via a namespace (and not attached):
#  [1] SummarizedExperiment_1.30.2 gtable_0.3.4                rjson_0.2.21                xfun_0.40                  
#  [5] Biobase_2.60.0              lattice_0.21-9              tzdb_0.4.0                  vctrs_0.6.4                
#  [9] tools_4.3.0                 bitops_1.0-7                generics_0.1.3              stats4_4.3.0               
# [13] parallel_4.3.0              fansi_1.0.5                 pkgconfig_2.0.3             Matrix_1.6-1.1             
# [17] S4Vectors_0.38.2            lifecycle_1.0.3             GenomeInfoDbData_1.2.10     compiler_4.3.0             
# [21] Rsamtools_2.16.0            Biostrings_2.68.1           munsell_0.5.0               codetools_0.2-19           
# [25] GenomeInfoDb_1.36.3         htmltools_0.5.6.1           RCurl_1.98-1.12             yaml_2.3.7                 
# [29] pillar_1.9.0                crayon_1.5.2                BiocParallel_1.34.2         DelayedArray_0.26.7        
# [33] abind_1.4-5                 tidyselect_1.2.0            digest_0.6.33               stringi_1.7.12             
# [37] restfulr_0.0.15             fastmap_1.1.1               grid_4.3.0                  colorspace_2.1-0           
# [41] cli_3.6.1                   magrittr_2.0.3              S4Arrays_1.0.6              XML_3.99-0.14              
# [45] utf8_1.2.3                  withr_2.5.1                 scales_1.2.1                timechange_0.2.0           
# [49] rmarkdown_2.25              XVector_0.40.0              matrixStats_1.0.0           gridExtra_2.3              
# [53] hms_1.1.3                   evaluate_0.22               knitr_1.44                  GenomicRanges_1.52.1       
# [57] IRanges_2.34.1              BiocIO_1.10.0               rtracklayer_1.60.1          rlang_1.1.1                
# [61] glue_1.6.2                  BiocManager_1.30.22         BiocGenerics_0.46.0         rstudioapi_0.15.0          
# [65] R6_2.5.1                    MatrixGenerics_1.12.3       GenomicAlignments_1.36.0    zlibbioc_1.46.0   