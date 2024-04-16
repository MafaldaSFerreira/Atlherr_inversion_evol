# Local PCA to genotype inversions

We run local PCAs within the inversion regions using lostruct.

First we generate bcf files that we can run in the program.

We use filtered genotypes maf 0.01 and we stick to 35 individuals from the Baltic Sea and Celtic Sea

~~~bash
#!/bin/bash -l
 
ml load bioinfo-tools bcftools/1.12

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/PCA/dataset_35_individuals"
input_VCF_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/filtered_vcf_miss20_maf0.01"

chromosomes=("chr6" "chr12" "chr17" "chr23")
coordinates=(6 12 17 23)
start=(20282765 15826318 23805445 14226443)
end=(26868581 27603093 29568511 19604273)

chrName=${chromosomes[$SLURM_ARRAY_TASK_ID]} 
c=${coordinates[$SLURM_ARRAY_TASK_ID]} 
s=${start[$SLURM_ARRAY_TASK_ID]}
e=${end[$SLURM_ARRAY_TASK_ID]}

cd ${WD}

### retrieve a vcf file for each inversion
bcftools view -S dataset_35_individuals.txt ${input_VCF_DIR}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${chrName}.miss20.maf0.01.vcf.gz | bcftools filter -e 'AC==0 || AC==AN' | bcftools view -m2 -M2 -v snps -O z -o herring_sentieon_91ind_190521.snps.inversion.${chrName}.d35.vcf.gz

tabix herring_sentieon_91ind_190521.snps.inversion.${chrName}.d35.vcf.gz

bcftools view -o ${chrName}.inv.d35.vcf herring_sentieon_91ind_190521.snps.inversion.${chrName}.d35.vcf.gz ${chrName}:${s}-${e}

bcftools convert -O b ${chrName}.inv.d35.vcf > ${chrName}.inv.d35.bcf

bcftools index ${chrName}.inv.d35.bcf
~~~

Run losttruct localy:

~~~bash
/Library/Frameworks/R.framework/Resources/Rscript ~/Documents/Postdoc/Repositories/local_pca/templated/run_lostruct.R -i chr12 -t snp -s 200 -I sample_info_all_SUK_all_Baltic.tsv -j 0007
/Library/Frameworks/R.framework/Resources/Rscript ~/Documents/Postdoc/Repositories/local_pca/templated/run_lostruct.R -i chr6 -t snp -s 200 -I sample_info_all_SUK_all_Baltic.tsv -j 0008
/Library/Frameworks/R.framework/Resources/Rscript ~/Documents/Postdoc/Repositories/local_pca/templated/run_lostruct.R -i chr17 -t snp -s 200 -I sample_info_all_SUK_all_Baltic.tsv -j 0009
/Library/Frameworks/R.framework/Resources/Rscript ~/Documents/Postdoc/Repositories/local_pca/templated/run_lostruct.R -i chr23 -t snp -s 200 -I sample_info_all_SUK_all_Baltic.tsv -j 00010
~~~

Check `localPCA.R` for the R code necessary to reconstruct the PCA scans in Figure ####

The input files and R session is available [link]
