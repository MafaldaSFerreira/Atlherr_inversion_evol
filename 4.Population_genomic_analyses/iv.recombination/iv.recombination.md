# Recombination R2

Start by creating a vcf file where we apply a more stringent filter of genotypes to save us some time in R2 calculations. We start from the set of genotypes with maf > 0.1 and miss < 20%. We filter genotypes with missing data higher than 10% and genotype quality lower than 30. 

### Generate vcf files

run_obtain_input_vcf_v2.sh
~~~bash
#!/bin/bash -l

ml load bioinfo-tools bcftools/1.12

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/recombination/repetition_02032023"
INPUT_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/filtered_vcf_miss20maf0.1"
OUTPUT_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/recombination/repetition_02032023/input_vcf_files"

# Chromosomes:
chromosomes=("chr6" "chr12" "chr17" "chr23" "chr1")

# obtain the file
ChrName="${chromosomes[$SLURM_ARRAY_TASK_ID]}"

bcftools filter -e 'F_MISSING > 0.1 && GQ < 30' ${INPUT_DIR}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.miss20.maf0.1.vcf.gz | bcftools view -m2 -M2 -v snps -O z -o ${OUTPUT_DIR}/herring_sentieon_91ind_190521.snps.miss90.GQ30.${ChrName}.vcf.gz
~~~

### Run vcftools

run_R2_snps.miss90.GQ30.full_chromosome.sh
~~~bash
#!/bin/bash -l
 
ml load bioinfo-tools vcftools/0.1.16

# Chromosomes:
chromosomes=("chr6" "chr12" "chr17" "chr23")
INPUT_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/recombination/repetition_02032023/input_vcf_files"
OUTPUT_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/recombination/repetition_02032023/R2_vcf.snps.miss90GQ30"
DATASETS="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/recombination/repetition_02032023/datasets"

# obtain the file
ChrName="${chromosomes[$SLURM_ARRAY_TASK_ID]}"

# LD Decay within Homozygote Classes:
vcftools --gzvcf ${INPUT_DIR}/herring_sentieon_91ind_190521.snps.miss90.GQ30.${ChrName}.vcf.gz --keep ${DATASETS}/homozygotes_${ChrName}_North.txt --geno-r2 --thin 10000 --out ${OUTPUT_DIR}/homozygotes_North_${ChrName}.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2

vcftools --gzvcf ${INPUT_DIR}/herring_sentieon_91ind_190521.snps.miss90.GQ30.${ChrName}.vcf.gz --keep ${DATASETS}/homozygotes_${ChrName}_South.txt --geno-r2 --thin 10000 --out ${OUTPUT_DIR}/homozygotes_South_${ChrName}.snps.miss90.GQ30.allchr.thin10kb.noLdwd.R2

# LD Decay within Homozygote Classes:
vcftools --gzvcf ${INPUT_DIR}/herring_sentieon_91ind_190521.snps.miss90.GQ30.${ChrName}.vcf.gz --keep ${DATASETS}/homozygotes_${ChrName}.txt --geno-r2 --thin 10000 --out ${OUTPUT_DIR}/only_homozygotes_${ChrName}.snps.miss90.GQ30.thin5kb.noLdwd.R2
~~~

### Plotting 

Check the R script `recombination_heatmaps.R` for the code to plot Supplementary Fig. 11
Check the R script `joint_R2_Fst_pi_plots.R` for the code to plot Figure 5 joint plots with Fst, Pi and R2 for each inversion.

The files to generate the plots and the .RData session can be found [link].

