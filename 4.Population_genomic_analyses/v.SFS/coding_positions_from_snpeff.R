#!/usr/bin/env Rscript
# This is an R script to extract non-coding positions from a SNPEff output
# ## Annotate SNPs:
# java -jar /proj/snic2020-2-19/private/herring/users/mafalda/software/snpEff/snpEff.jar Ch_v2.0.2.105 atlantic_herring_61i.${prefix}.inv.${ChrName}.vcf.gz > atlantic_herring_61i.${prefix}.inv.${ChrName}.ann.vcf
# Extract them using snpsift
# bcftools view -v snps -H atlantic_herring_61i.miss20.inv.chr12.ann.vcf.gz | /proj/snic2020-2-19/private/herring/users/mafalda/software/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /proj/snic2020-2-19/private/herring/users/mafalda/software/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT ANN[*].GENEID ANN[*].FEATUREID ANN[*].IMPACT ANN[*].EFFECT > annotations/atlantic_herring_61i.miss20.inv.chr12.ann.ANNO.txt
# For more details, check "9.Mutational_load.md" on date: ### 2022-10-19 

# Mafalda Ferreira, 2022-10-22
# usage: Rscript coding_positions_from_snpeff.R input_file chromosome output_nonsyn_file output_syn_file

# Libraries
suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library("argparse"))

# Read arguments and input file name
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
chromosome <- as.numeric(args[2])
output_nonsyn_file <- args[3]
output_syn_file <- args[4]

# parser <- ArgumentParser(description='Extract non-synonymous and synonymous positions from SNPEff annotations.')
# parser$add_argument(arg="--input", help="input file")
# parser$add_argument(arg="--chromosome", help="chromosome being processed", type="integer")
# parser$add_argument(arg="--nonsyn", help="name of output of nonsynonymous positions", default="nonsynonymous_positions.txt")
# parser$add_argument(arg="--syn", help="name of output of synonymous positions", default="synonymous_positions.txt")
# 
# args <- parser$parse_args()

# Set working directory
setwd("./")
WD=getwd()

# Read files
annotation<-read.table(input_file, skip=1)

colnames(annotation)<-c("Chr", "Pos", "Ref", "Alt", "GeneID", "TranscriptID","Impact", "Effect")

# This comes from Mats code, but in this case we just annotation non-synonymous and synonymous mutations
annotation[,"NonSyn"]<-grepl("missense_variant", annotation$Effect)
annotation[,"Syn"] <- grepl("synonymous_variant", annotation$Effect)

# What are all coding positions
annotation_coding <- annotation %>% filter(NonSyn == T | Syn==T)
print(paste("Number of total nonsynonymous + synonymous positions: ", length(unique(annotation_coding$Pos))))

# Extract non_syn positions
annotation_nonsyn_pos <- annotation %>% filter(NonSyn == T ) %>% distinct(Pos)
# How many non-synonymous positions?
print(paste0("Number of nonsynonymous positions in this chromosomes: ", length(unique(annotation_nonsyn_pos$Pos))))


# Now, let's extract the synonymous positions, making sure they ARE Not already in the nonsyn_positions. This will avoid duplicated positions
annotation_syn_pos <- annotation_coding %>% filter(Syn == T & !Pos %in% annotation_nonsyn_pos$Pos) %>% distinct(Pos)
# How many synonymous positions?
print(paste0("Number of synonymous positions in this chromosomes: ", length(unique(annotation_syn_pos$Pos))))

# Is the length of the coding positions the same?
print(paste("Is the length of the coding positions the same?",length(unique(annotation_syn_pos$Pos))+length(unique(annotation_nonsyn_pos$Pos)) == length(unique(annotation_coding$Pos))))

# Let's write the outputs
write.table(data.frame(CHR=chromosome, POS=annotation_nonsyn_pos$Pos), col.names = F, row.names = F, quote = F, sep="\t", file = output_nonsyn_file)

write.table(data.frame(CHR=chromosome, POS=annotation_syn_pos$Pos), col.names = F, row.names = F, quote = F, sep="\t", file = output_syn_file)

