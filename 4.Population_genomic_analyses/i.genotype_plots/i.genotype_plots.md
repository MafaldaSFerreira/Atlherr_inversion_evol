# Genotype plots

Extract genotypes for highly differentiated SNPs from the vcf file:

First, split the file 
~~~bash
cut -f1,2 -d' ' NSvNE_logp15.txt | awk '{print $1 "\t" $2}' > NSvNE_logp15.SNPs.txt

grep "^chr6" NSvNE_logp15.SNPs.txt > NSvNE_logp15.SNPs.chr6.txt
grep "^chr12" NSvNE_logp15.SNPs.txt > NSvNE_logp15.SNPs.chr12.txt
grep "^chr17" NSvNE_logp15.SNPs.txt > NSvNE_logp15.SNPs.chr17.txt
grep "^chr23" NSvNE_logp15.SNPs.txt > NSvNE_logp15.SNPs.chr23.txt
~~~

Extract the genotypes

~~~bash
#!/bin/bash -l

ml load bcftools/1.12

chromosomes=("chr6" "chr12" "chr17" "chr23")

# obtain the Chr
ChrName="${chromosomes[$SLURM_ARRAY_TASK_ID]}"

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/Genotypes"
INPUT_VCFS_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/filtered_vcf_miss20maf0.1"

cd ${WD}

#filter vcf file
bcftools view --regions-file NSvNE_logp15.SNPs.${ChrName}.txt -Oz -o herring_sentieon_91ind_190521.${ChrName}.NSvNE_logp15.vcf.gz ${INVPUT_VCFS_DIR}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.miss20.maf0.1.vcf.gz  

bcftools view -m2 -M2 -v snps herring_sentieon_91ind_190521.${ChrName}.NSvNE_logp15.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' > herring_sentieon_91ind_190521.${ChrName}.NSvNE_logp15.nomono.nomulti.genotypes

sed 's/|/\t/g' herring_sentieon_91ind_190521.${ChrName}.NSvNE_logp15.nomono.nomulti.genotypes > herring_sentieon_91ind_190521.${ChrName}.NSvNE_logp15.nomono.nomulti.sep.genotypes
~~~


Check `genotype_heatmaps.R` to see code to plot the heatmap in Figure ####

