# Site Frequency Spectrum for non-coding and coding sites

### SNPEff annotation

First, we annotate a VCF file with all 61 Atlantic herring individuals using SNPEff

~~~bash
## Load required modules
ml load bioinfo-tools bcftools/1.12
ml load java/OracleJDK_11.0.9

## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

# What is the prefix to name the vcf file:
prefix="miss20"

# Working directory
INPUT_VCF="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/filtered_vcf_miss20_maf0.01"
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/SFS"
software_WD="/proj/snic2020-2-19/private/herring/users/mafalda/software"
vcf_files_DIR=${WD}"/vcf_files_61i"
annotations_DIR=${WD}"/annotations"

## Retrieve Atlantic herring samples:
#echo "Generating vcf files..."

bcftools view -S ${WD}/sample_files/samples_atlantic_herring_61_individuals.txt -R /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/MKT_Sprat/inputs/alignments_24022023/${ChrName}.geno.files.24022023.chr.bed ${INPUT_VCF}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.miss20.maf0.01.vcf.gz -O z -o ${vcf_files_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.vcf.gz

bcftools index ${vcf_files_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.vcf.gz

## Annotate SNPs:
echo "Annotating vcf files..."

java -jar ${software_WD}/snpEff/snpEff.jar Ch_v2.0.2.105 ${vcf_files_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.vcf.gz > ${vcf_files_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.ann.vcf

bgzip ${vcf_files_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.ann.vcf
bcftools index ${vcf_files_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.ann.vcf.gz

## Extract SNPEff annotation:
echo "Extracting coding positions..."

bcftools view -v snps -H ${vcf_files_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.ann.vcf.gz | ${software_WD}/snpEff/scripts/vcfEffOnePerLine.pl | java -jar ${software_WD}/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT ANN[*].GENEID ANN[*].FEATUREID ANN[*].IMPACT ANN[*].EFFECT > ${annotations_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.ann.ANNO.txt

#echo "Let's go to R..."
# Load R packages here, because it will re-load a version of java incompatible with SNPEff
ml load R_packages

Rscript ${software_WD}/coding_positions_from_snpeff.R ${annotations_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.ann.ANNO.txt ${ChrName/chr/} ${annotations_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.ann.ANNO.nonsyn.txt ${annotations_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.ann.ANNO.syn.txt

rm ${vcf_files_DIR}/atlantic_herring_61i.${prefix}.inv.${ChrName}.vcf.gz
~~~

### Add ancestral state to the vcf file

Here, we use a vcf file not filtered for minor allele frequency, since the SFS needs to include all sites.

First, we merge our files with the sprat vcf file. See how I generated the Sprat vcf file in `2.Mapping_sprat_assembly`

~~~bash
## Load required modules
ml load bioinfo-tools bcftools/1.12

## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

# What is the prefix to name the vcf file:
prefix="miss20"

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/SFS/repetition_02032023"
SPRAT_VCF="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/MKT_Sprat/vcf_files_alignments_0.25p_24022023/"
OUTPUT_VCF="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/SFS/repetition_02032023/intermediate_vcf"

## Filter samples 

cd ${WD}

bcftools view -e 'F_MISSING > 0.2' -O z -o ${OUTPUT_VCF}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${prefix}.${ChrName}.vcf.gz herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.vcf.gz

bcftools index ${OUTPUT_VCF}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${prefix}.${ChrName}.vcf.gz

bcftools merge ${OUTPUT_VCF}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${prefix}.${ChrName}.vcf.gz ${SPRAT_VCF}/Ch_EuSprat.${ChrName}.sorted.norm.vcf.gz | bcftools annotate --rename-chrs /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/MKT_Sprat/map.file -O z -o ${OUTPUT_VCF}/sprat_herring_91i.${prefix}.${ChrName}.renamed.vcf.gz

tabix ${OUTPUT_VCF}/sprat_herring_91i.${prefix}.${ChrName}.renamed.vcf.gz
~~~

Inject the ancestral state in the INFO field of the vcf file

~~~bash
# load modules
ml load bioinfo-tools vcftools bcftools

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/SFS/repetition_02032023/intermediate_vcf"

cd ${WD}

# determine Chr
ChrName=chr${SLURM_ARRAY_TASK_ID}

bcftools view -s EuSprat sprat_herring_91i.miss20.${ChrName}.renamed.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%IUPACGT\t]\n' > ancestral_allele_${ChrName}.txt

sed 's/\.\/\./\./' ancestral_allele_${ChrName}.txt > ancestral_allele_miss_${ChrName}.txt
    
# create a header file with the information for the header: ancestral_allele_header.hdr
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">

# index specifying that the start and end are the second column
bgzip ancestral_allele_miss_${ChrName}.txt
tabix -s1 -b2 -e2 -f ancestral_allele_miss_${ChrName}.txt.gz

bcftools annotate -O z -o sprat_atlantic_herring_91i.${ChrName}.AA.vcf.gz -c CHROM,POS,+AA -a ancestral_allele_miss_${ChrName}.txt.gz -h ancestral_allele_header.hdr sprat_herring_91i.miss20.${ChrName}.renamed.vcf.gz
~~~

### Calculate derived allele frequencies

Calculate derived allele frequencies. I generated files where I state which individuals are homozygotes for each inversion. 

~~~bash
# NOTE: This script assumes the input vcf files are already generated. Check 2022-09-07, step 2, to see how we have generated these files, and scripts run_generate_vcf.sh and run_ancestral.sh on the working directory scripts folder.

## Load required modules
ml load bioinfo-tools vcftools/0.1.16

# Chromosomes
chromosomes=("chr6" "chr12" "chr17" "chr23")
ChrName="${chromosomes[$SLURM_ARRAY_TASK_ID]}"

# Working directory
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/SFS/repetition_02032023"
vcf_files_DIR=${WD}"/intermediate_vcf"
annotations_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/SFS/annotations"
freq_results=${WD}"/output_derived_frequencies"
samples_folder=${WD}"/datasets"


# Celtic Sea individuals
vcftools --keep ${samples_folder}/homozygotes_${ChrName}_South.txt --max-alleles 2 --min-alleles 2 --positions ${annotations_DIR}/atlantic_herring_61i.miss20.inv.${ChrName}.ann.ANNO.nonsyn.txt --gzvcf ${vcf_files_DIR}/sprat_atlantic_herring_91i.${ChrName}.AA.vcf.gz --freq --derived --out ${freq_results}/frequencies_South_${ChrName}_nonsyn
vcftools --keep ${samples_folder}/homozygotes_${ChrName}_South.txt --max-alleles 2 --min-alleles 2 --positions ${annotations_DIR}/atlantic_herring_61i.miss20.inv.${ChrName}.ann.ANNO.syn.txt --gzvcf ${vcf_files_DIR}/sprat_atlantic_herring_91i.${ChrName}.AA.vcf.gz --freq --derived --out ${freq_results}/frequencies_South_${ChrName}_syn

sed 's/:/\t/g' ${freq_results}/frequencies_South_${ChrName}_nonsyn.frq > ${freq_results}/frequencies_South_${ChrName}_nonsyn.mod.frq
sed 's/:/\t/g' ${freq_results}/frequencies_South_${ChrName}_syn.frq > ${freq_results}/frequencies_South_${ChrName}_syn.mod.frq

# Baltic Sea individuals:
vcftools --keep ${samples_folder}/homozygotes_${ChrName}_North.txt --max-alleles 2 --min-alleles 2 --positions ${annotations_DIR}/atlantic_herring_61i.miss20.inv.${ChrName}.ann.ANNO.nonsyn.txt --gzvcf ${vcf_files_DIR}/sprat_atlantic_herring_91i.${ChrName}.AA.vcf.gz --freq --derived --out ${freq_results}/frequencies_North_${ChrName}_nonsyn
vcftools --keep ${samples_folder}/homozygotes_${ChrName}_North.txt --max-alleles 2 --min-alleles 2 --positions ${annotations_DIR}/atlantic_herring_61i.miss20.inv.${ChrName}.ann.ANNO.syn.txt --gzvcf ${vcf_files_DIR}/sprat_atlantic_herring_91i.${ChrName}.AA.vcf.gz --freq --derived --out ${freq_results}/frequencies_North_${ChrName}_syn

sed 's/:/\t/g' ${freq_results}/frequencies_North_${ChrName}_nonsyn.frq > ${freq_results}/frequencies_North_${ChrName}_nonsyn.mod.frq
sed 's/:/\t/g' ${freq_results}/frequencies_North_${ChrName}_syn.frq > ${freq_results}/frequencies_North_${ChrName}_syn.mod.frq
~~~

Calculate derived allele frequencies for the entire genome, to compare with the inversion regions:

~~~bash
## Load required modules
ml load bioinfo-tools vcftools/0.1.16

# Chromosomes
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Working directory
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/SFS/repetition_02032023"
vcf_files_DIR=${WD}"/intermediate_vcf"
annotations_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Mutational_load_SFS/SFS/annotations"
freq_results=${WD}"/output_derived_frequencies"

# All individuals:
vcftools --max-alleles 2 --min-alleles 2 --positions ${annotations_DIR}/atlantic_herring_61i.miss20.inv.${ChrName}.ann.ANNO.nonsyn.txt --gzvcf ${vcf_files_DIR}/sprat_atlantic_herring_91i.${ChrName}.AA.vcf.gz --freq --derived --out ${freq_results}/frequencies_${ChrName}_nonsyn
vcftools  --max-alleles 2 --min-alleles 2 --positions ${annotations_DIR}/atlantic_herring_61i.miss20.inv.${ChrName}.ann.ANNO.syn.txt --gzvcf ${vcf_files_DIR}/sprat_atlantic_herring_91i.${ChrName}.AA.vcf.gz --freq --derived --out ${freq_results}/frequencies_${ChrName}_syn

sed 's/:/\t/g' ${freq_results}/frequencies_${ChrName}_nonsyn.frq > ${freq_results}/frequencies_North_${ChrName}_nonsyn.mod.frq
sed 's/:/\t/g' ${freq_results}/frequencies_${ChrName}_syn.frq > ${freq_results}/frequencies_North_${ChrName}_syn.mod.frq
~~~

## Plotting

Use code in `plotting_SFS.R` do generate plots in Figure 6 and Supplementary Figure 13.

The input files used in the R code and output figures can be found in [link]