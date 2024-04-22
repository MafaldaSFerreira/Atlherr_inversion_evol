# Ancestral state from sprat genome

We start with the output of Satsuma and use R code to extract regions in the sprat assembly that map to genes in the Atlantic herring genome.

~~~R
library(data.table)
library(biomaRt)
library(GenomicRanges)
library(tidyverse)
require(Biostrings)
require(GenomicRanges)
require(ape)

# Load annotations from Ensembl ####

ensembl <- useEnsembl(biomart="ensembl",dataset="charengus_gene_ensembl", mirror="www")

attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

genes <- getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position"), mart = ensembl)

# Load satsuma output and the reference genomes ####

IPA_pri_satsuma <- read.table("/Users/maffe217/Dropbox/Mac/Documents/Postdoc/Herring/Inversion_Project/Satsuma_Align_Mats/hifiasm_assemblies/satsuma_summary.chained.out_sorted_sprat_hap1", stringsAsFactors=F, sep = "\t", comment.char = "")
IPA_pri <- readDNAStringSet("/Users/maffe217/Dropbox/Mac/Documents/Postdoc/Herring/Inversion_Project/Satsuma_Align_Mats/hifiasm_assemblies/sprat.hap1.fa")
Herring <- readDNAStringSet("/Users/maffe217/Dropbox/Mac/Documents/Postdoc/Herring/Assembly/Ch_v2.0.2.fa")

# Create an output matrix, where we will obtain the section of the sprat region that is aligning to a gene. 

output<-data.frame(matrix(data=NA,nrow=1,ncol=8))

colnames(output)<-c("gene","chromosome","start","end","sprat_contig","sprat_start","sprat_end","sprat_strand")

for(i in 1:nrow(genes)){

  gene <- genes[i,]
  chr=paste0("chr",gene[,2])
  start=gene[,3] 
  end=gene[,4] 
  
  # here we extract the section of the alignment where this gene is contained.
  
  IPA_pri_satsuma[IPA_pri_satsuma$V4 == chr & IPA_pri_satsuma$V6 >= start & IPA_pri_satsuma$V5 <= end,]->tmp
  
  # extract the columns we want and name them
  tmp2<-tmp[,c(1,2,3,8)]
  colnames(tmp2)<-c("sprat_contig","sprat_start","sprat_end","sprat_strand")
  
  # Add herring gene information:
  tmp2$gene<-base::rep(gene[,1],nrow(tmp2))
  tmp2$chromosome<-base::rep(chr,nrow(tmp2))
  tmp2$start<-base::rep(start,nrow(tmp2))
  tmp2$end<-base::rep(end,nrow(tmp2))
  
  # Add these to output
  rbind(output,tmp2)->output
  
}

# Write output to a file
output<-na.omit(output)
write.table(output, file="gene_mappings_herring_vs_sprat_14032022.txt", col.names =T, row.names = F, quote=F, sep="\t")

# It seems stasuma coordinates are zero based... so let's try just adding 1
output$sprat_start<-output$sprat_start+1
output$sprat_end<-output$sprat_start+1

# Summarize:
output %>%
  group_by(gene,sprat_contig)%>%
  summarise(start=min(sprat_start),
            end=max(sprat_end),
            strand=unique(sprat_strand)) %>%
  mutate(alig_len= (end-start)) -> output_by_gene_with_length

# choose the longest alignment if more than 1 contigs map to one gene
# top_n will chose the longest. but if the lenght of the mapping is the same
# for two contigs, top_n will keep all of them
# distinct will select the first entry of those by gene.
output_by_gene_with_length %>% 
  group_by(gene) %>%
  top_n(1, alig_len) %>%
  distinct(gene, .keep_all=T) -> output_by_gene_with_length_unique

nrow(output_by_gene_with_length_unique)
#23911 #gene_mappings_herring_vs_sprat_14032022.txt

merge(genes,output_by_gene_with_length_unique, by.x="ensembl_gene_id",by.y="gene") -> mappings_1_vs_1

length(unique(mappings_1_vs_1$ensembl_gene_id)) 
# [1] 23911 #gene_mappings_herring_vs_sprat_14032022.txt

mappings_1_vs_1$gene_length<-mappings_1_vs_1$end_position-mappings_1_vs_1$start_position

mappings_1_vs_1$prop_align<-mappings_1_vs_1$alig_len/mappings_1_vs_1$gene_length

plot(mappings_1_vs_1$gene_length,mappings_1_vs_1$prop_align)

ggplot(mappings_1_vs_1) +
  geom_density(aes(x=prop_align))+
  xlim(0,5)

# ok, let's filter these on some random length.
subset(mappings_1_vs_1,mappings_1_vs_1$prop_align>=0.25)->mappings_1_vs_1_0.25p

# how many genes we have left?
length(unique(mappings_1_vs_1_0.25p$ensembl_gene_id))
# [1] 22635 #gene_mappings_herring_vs_sprat_14032022.txt


# Using the coordinates of the alignment, extract a fasta file containing the herring gene and the respective sprat region mapping to it. 
mappings_1_vs_1_0.25p$chromosome_name<-paste0("chr",mappings_1_vs_1_0.25p$chromosome_name)

mappings_1_vs_1_0.25p[mappings_1_vs_1_0.25p$start==0,]

mappings_1_vs_1_0.25p$sprat_chromosome<-str_split_fixed(mappings_1_vs_1_0.25p$sprat_contig, "_", 2)[,1]

for(i in 1:nrow(mappings_1_vs_1_0.25p)){
  
  g=mappings_1_vs_1_0.25p[i,"ensembl_gene_id"]
  
  herring_GR = GRanges(paste0(mappings_1_vs_1_0.25p[i,2],":",mappings_1_vs_1_0.25p[i,3],"-",mappings_1_vs_1_0.25p[i,4]))
  sprat_pri_GR = GRanges(paste0(mappings_1_vs_1_0.25p[i,"sprat_chromosome"],":",mappings_1_vs_1_0.25p[i,6],"-",mappings_1_vs_1_0.25p[i,7]))
  
  output_name=paste0("./input_fasta_0.25p/",g,"_sprat_0.25p.fa")
  
  writeXStringSet(x = Herring[herring_GR], file = output_name, append=T)
  
  
  if(mappings_1_vs_1_0.25p[i,8]=="-"){

    writeXStringSet(x = reverseComplement(IPA_pri[sprat_pri_GR]), file = output_name, append = T)
    
  }else{
    
    writeXStringSet(x = IPA_pri[sprat_pri_GR], file = output_name, append=T)
  }

}

# Create info files that contain the ENSEMBL gene ID, chromosome, start and end for each gene.
for(i in 1:nrow(mappings_1_vs_1_0.25p)){

  tmp<-mappings_1_vs_1_0.25p[i,c(1,2,3,4)]
  # fix chromosome names
  tmp$chromosome_name<-str_remove(tmp$chromosome_name, "chr")
  
  write.table(tmp,file=paste0("./input_info_fasta_0.25p/",tmp$ensembl_gene_id,".info"),sep="\t",row.names=F,col.names=F,quote=F)
  
}
~~~

### Realign with Mafft
Now, we align these genes using mafft and AMAS to obtain a summary of the alignments` quality.

~~~bash
ml load bioinfo-tools MAFFT/7.407

parallel -j 20 'i={1}; mafft --quiet --auto $i > ${i/.fa/.algn.fa}' ::: $(ls *0.25p.fa)

parallel -j 18 'i={1}; python /proj/snic2020-2-19/private/herring/users/mafalda/software/AMAS/amas/AMAS.py summary -i $i -f fasta -d dna -o ${i/.fa/.summary.txt}' ::: $(ls *0.25p.algn.v2.fa)
~~~

Concate the summaries, to obtain an table with all information and select the alingments where missing data is bellow 20 and the proportion of variable sites is bellow 0.2.

~~~bash
ml load R_packages

wg_tables.R summary.txt algn.v2.summary.txt

cut -f1,6,8,10 maff_alignments_version2.summary.txt | awk '{if($2 <= 20 && $3 <=0.2) print $0}' | cut -f1 > maff_alignments_version2.Miss20VariableSites20.txt
~~~

### Convert fasta to geno
We convert the alingments to .geno files, a format compatible with [`Genomics General`](https://github.com/simonhmartin/genomics_general/blob/master/VCF_processing/genoToVCF.py) script genoToVCF.py, which we will use to convert the geno files into vcf files.

~~~bash
parallel -j 20 'i={1}; python3 /proj/snic2020-2-19/private/herring/users/mafalda/software/ancestral_vcf.py $i ${i/_sprat_0.25p.algn.fa/.info} ${i/_sprat_0.25p.algn.fa/.geno}' ::: $(ls *0.25p.algn.fa)
~~~

### Convert geno into vcfs

In R, I generated the files ".geno.files.14032022.txt" for each chromosome, which I will use to concatenate the genos files using genomics generals. These files are sorted by start position, so that we have a sorted geno file per chromosome when we merge.

~~~R
output<-read.table("gene_mappings_herring_vs_sprat_14032022.txt",header=T,sep="\t")

head(output)
length(unique(output$gene))
23911 # number of genes covered

# It seems stasuma coordinates are zero based... so let's try just adding 1
output$sprat_start<-output$sprat_start+1
output$sprat_end<-output$sprat_start+1

# Summarize:
output %>%
  group_by(gene,sprat_contig)%>%
  summarise(start=min(sprat_start),
            end=max(sprat_end),
            strand=unique(sprat_strand)) %>%
  mutate(alig_len= (end-start)) -> output_by_gene_with_length

# chose the longest alignment if more than 1 contigs map to one gene
# top_n will chose the longest. but if the lenght of the mapping is the same
# for two contigs, top_n will keep all of them
# distinct will select the first entry of those by gene.
output_by_gene_with_length %>% 
  group_by(gene) %>%
  top_n(1, alig_len) %>%
  distinct(gene, .keep_all=T) -> output_by_gene_with_length_unique

nrow(output_by_gene_with_length_unique)
23911 #gene_mappings_herring_vs_sprat_14032022.txt

# so, now I think I should probably remove hits that are too short.

merge(genes,output_by_gene_with_length_unique, by.x="ensembl_gene_id",by.y="gene") -> mappings_1_vs_1

length(unique(mappings_1_vs_1$ensembl_gene_id)) 
# [1] 23911 #gene_mappings_herring_vs_sprat_14032022.txt

mappings_1_vs_1$gene_length<-mappings_1_vs_1$end_position-mappings_1_vs_1$start_position

mappings_1_vs_1$prop_align<-mappings_1_vs_1$alig_len/mappings_1_vs_1$gene_length

plot(mappings_1_vs_1$gene_length,mappings_1_vs_1$prop_align)

ggplot(mappings_1_vs_1) +
  geom_density(aes(x=prop_align))+
  xlim(0,5)

# ok, let's filter these on some random length.
subset(mappings_1_vs_1,mappings_1_vs_1$prop_align>=0.25)->mappings_1_vs_1_0.25p

# how many genes we have left?
length(unique(mappings_1_vs_1_0.25p$ensembl_gene_id))
# [1] 22635 #gene_mappings_herring_vs_sprat_14032022.txt


# ok, and now essentially I just need to generate
mappings_1_vs_1_0.25p$chromosome_name<-paste0("chr",mappings_1_vs_1_0.25p$chromosome_name)

mappings_1_vs_1_0.25p[mappings_1_vs_1_0.25p$start==0,]

#Now I want to filter these genes on variale sites and missing data:

filter<-read.table("/Users/maffe217/Documents/Postdoc/Project_Herring/Inversion_Project/Satsuma_Align_Mats/hifiasm_assemblies/maff_alignments_version2.Miss20VariableSites20.txt")

genes_to_keep<-str_remove(filter[,1],"_sprat_0.25p.algn.fa")

### output a merged file.

mappings_1_vs_1_0.25p %>%
  filter(ensembl_gene_id %in% genes_to_keep) %>%
  group_by(chromosome_name, ensembl_gene_id) %>%
  summarise(start=min(start_position), end=max(end_position)) -> df

df$chromosome_name<-paste0("chr", df$chromosome_name)

#nrow(df)
#[1] 15471

for(i in 1:length(unique(df$chromosome_name))){
  
  chr=unique(df$chromosome_name)[i]
  
  tmp<-subset(df, df$chromosome_name == chr)
  # sort by start position!
  tmp<-tmp[order(tmp$start),]
  
  tmp$file<-paste0(tmp$ensembl_gene_id,".geno")
  
  write.table(tmp,file=paste0(chr,".geno.files.24022023.txt"),quote = F, col.names = F, row.names = F, sep="\t")
}
~~~

And now, we use these geno files and the code bellow to create a vcf file for each chromosome with the sprat genotypes.

~~~bash
#!/bin/bash -l
 
# Load packages
ml load R_packages
ml load gnuparallel 
ml load python/3.8.7
ml load bioinfo-tools bcftools/1.12

# Chromosome name
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Merge the geno files
# Start process inside folder alignments_0.25p_010322/
cd /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/MKT_Sprat/alignments_0.25p_010322

Rscript /proj/snic2020-2-19/private/herring/users/mafalda/software/merge_gene_genos.R /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/MKT_Sprat/inputs/alignments_24022023/"${ChrName}".geno.files.24022023.txt vcf_files_alignments_0.25p_24022023

# Move to the target folder
cd /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/MKT_Sprat/vcf_files_alignments_0.25p_24022023

# Convert the calls into a pseudo-phased format. Since the assembly is haploid, the calls aren't trully phased, and we will only have homozygotes
python3 /proj/snic2020-2-19/private/herring/users/mafalda/software/genomics_generals_Oct21/filterGenotypes.py -if diplo -of phased -i "${ChrName}".genes.geno -o "${ChrName}".genes.phased.geno -t 1

# Convert geno file to vcf file
python3 /proj/snic2020-2-19/private/herring/users/mafalda/software/genomics_generals_Oct21/genoToVCF.py -g "${ChrName}".genes.phased.geno -f phased -o Ch_EuSprat."${ChrName}".vcf -r /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/MKT_Sprat/reference/Ch_v2.0.2."${ChrName}".fasta

# Zip and index the vcf
bgzip Ch_EuSprat."${ChrName}".vcf

# Sort the vcf file
bcftools sort -O z -o Ch_EuSprat."${ChrName}".sorted.vcf.gz Ch_EuSprat."${ChrName}".vcf.gz

tabix Ch_EuSprat."${ChrName}".sorted.vcf.gz

# Remove duplicated positions
bcftools norm -O z --rm-dup all -o Ch_EuSprat."${ChrName}".sorted.norm.vcf.gz Ch_EuSprat."${ChrName}".sorted.vcf.gz

# Index vcf
tabix Ch_EuSprat."${ChrName}".sorted.norm.vcf.gz

# Merge file with the rest of the individuals
bcftools merge -O z -o herring_sprat."${ChrName}".vcf.gz /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv."${ChrName}".vcf.gz Ch_EuSprat."${ChrName}".sorted.norm.vcf.gz

# Rename chromosome names
bcftools annotate --rename-chrs /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/MKT_Sprat/map.file -O z -o herring_sprat."${ChrName}".renamed.vcf.gz herring_sprat."${ChrName}".vcf.gz

# Index
tabix herring_sprat."${ChrName}".renamed.vcf.gz

# Remove intermediate files
#rm "${ChrName}".genes.phased.geno
#rm Ch_EuSprat."${ChrName}".vcf.gz
#rm Ch_EuSprat."${ChrName}".sorted.vcf.gz
#rm herring_sprat."${ChrName}".vcf.gz

~~~

### Vcf to Consensus Fasta in the coordinates of the Atlantic herring genome

After getting a vcf file, we can now insert the sprat variants into the Atlantic herring genome and generate a whole-genome consensus fasta

~~~bash
#!/bin/bash -l
 
## Load required modules
ml load bioinfo-tools bcftools/1.12 BEDTools/2.29.2 samtools/1.12

## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/Sprat"
INPUT_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/MKT_Sprat/vcf_files_alignments_0.25p_24022023"

cd ${WD}

bcftools view -s EuSprat ${INPUT_DIR}/Ch_EuSprat.${ChrName}.sorted.norm.vcf.gz | bcftools view -H -g ^miss | awk '{OFS="\t"}{print $1,$2-1,$2}' > ${ChrName}.sprat.called.bed

grep -w ${ChrName} /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/Ch_v2.0.2.chromsizes > ${ChrName}.chromsize

bedtools complement -i ${ChrName}.sprat.called.bed -g ${ChrName}.chromsize > ${ChrName}.sprat.nocall.bed

samtools faidx /crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta ${ChrName} | bcftools consensus ${INPUT_DIR}/Ch_EuSprat.${ChrName}.sorted.norm.vcf.gz --sample EuSprat > EuSprat.${ChrName}.fasta

bedtools maskfasta -fi EuSprat.${ChrName}.fasta -bed ${ChrName}.sprat.nocall.bed -fo EuSprat.${ChrName}.consensus.fa

rm ${ChrName}.sprat.nocall.bed
rm ${ChrName}.sprat.called.bed
rm EuSprat.${ChrName}.fasta
~~~

These files are used for dNdS calculations, SFS and phylogenetic analyses.