# Pixy


### Divergence statistics between homozygotes for the Northern and Southern allele

I ran pixy for individuals ascertained as homozygote for the Northern or Southern allele to infer dxy and fst. Since the homozygotes vary between inversions, I have four different slurm scripts that use a different population file specific for each inversion. Check code bellow:

run_chr12_v02.sh
~~~bash
#!/bin/bash -l
 
## Load required modules
ml load bioinfo-tools pixy/1.2.5.beta1
## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Read arguments from the command line
DIR1=$ARG1
DIR2=$ARG2
FILTER=$ARG3

# Define directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}
PopDIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}"/population_files"
input_vcf_dir="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/"${DIR2}

cd ${WD}

echo ${input_vcf_dir}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.${ARG3}.vcf.gz

pixy --stats pi fst dxy --populations ${PopDIR}/chr12_population_file.txt --vcf ${input_vcf_dir}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.${ARG3}.vcf.gz --window_size 20000 --n_cores 2 --output_folder chr12_inversion_popgen_stats --output_prefix chr12_dataset42i_PacifHerr.${ChrName}.${ARG3}.20kb.popgenpixy.out
~~~

run_chr17_v02.sh
~~~bash
#!/bin/bash -l
 
## Load required modules
ml load bioinfo-tools pixy/1.2.5.beta1
## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Read arguments from the command line
DIR1=$ARG1
DIR2=$ARG2
FILTER=$ARG3

# Define directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}
PopDIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}"/population_files"
input_vcf_dir="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/"${DIR2}

cd ${WD}

pixy --stats pi fst dxy --populations ${PopDIR}/chr17_population_file.txt --vcf ${input_vcf_dir}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.${ARG3}.vcf.gz --window_size 20000 --n_cores 2 --output_folder chr17_inversion_popgen_stats --output_prefix chr17_dataset42i_PacifHerr.${ChrName}.${ARG3}.20kb.popgenpixy.out
~~~

run_chr6_v02.sh
~~~bash
#!/bin/bash -l

## Load required modules
ml load bioinfo-tools pixy/1.2.5.beta1
## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Read arguments from the command line
DIR1=$ARG1
DIR2=$ARG2
FILTER=$ARG3

# Define directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}
PopDIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}"/population_files"
input_vcf_dir="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/"${DIR2}

cd ${WD}

pixy --stats pi fst dxy --populations ${PopDIR}/chr6_population_file.txt --vcf ${input_vcf_dir}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.${ARG3}.vcf.gz --window_size 20000 --n_cores 2 --output_folder chr6_inversion_popgen_stats --output_prefix chr6_dataset42i_PacifHerr.${ChrName}.${ARG3}.20kb.popgenpixy.out
~~~

run_chr23_v02.sh
~~~bash
#!/bin/bash -l
 
## Load required modules
ml load bioinfo-tools pixy/1.2.5.beta1
## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Read arguments from the command line
DIR1=$ARG1
DIR2=$ARG2
FILTER=$ARG3

# Define directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}
PopDIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}"/population_files"
input_vcf_dir="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/"${DIR2}

cd ${WD}

pixy --stats pi fst dxy --populations ${PopDIR}/chr23_population_file.txt --vcf ${input_vcf_dir}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.${ARG3}.vcf.gz --window_size 20000 --n_cores 2 --output_folder chr23_inversion_popgen_stats --output_prefix chr23_dataset42i_PacifHerr.${ChrName}.${ARG3}.20kb.popgenpixy.out
~~~

### Divergence statistics between and Pacific herring

I also calculate statistics between all Atlantic herring (N=61) and all Pacific herring (N=30) in the dataset with the following script:
~~~bash
#!/bin/bash -l
 
## Load required modules
ml load bioinfo-tools pixy/1.2.5.beta1
## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Read arguments from the command line
DIR1=$ARG1
DIR2=$ARG2
FILTER=$ARG3

# Define directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}
PopDIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}"/population_files"
input_vcf_dir="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/"${DIR2}

# What are the arguments?
echo ${WD}
echo ${PopDIR}
echo ${input_vcf_dir}

cd ${WD}

pixy --stats pi fst dxy --populations ${PopDIR}/atlantic_vs_pacific_all_populations_file.txt --vcf ${input_vcf_dir}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.${ARG3}.vcf.gz --window_size 20000 --n_cores 4 --output_folder atlantic_vs_pacific_popgen_stats --output_prefix atlantic_vs_pacific.${ChrName}.20kb.popgenpixy.out
~~~

These slurm scripts were ran as follows:

~~~bash
sbatch --export=ALL,ARG1="repetition_miss20maf0.01",ARG2="filtered_vcf_miss20_maf0.01",ARG3="miss20.maf0.01" run_atl_vs_pac_v02.sh
sbatch --export=ALL,ARG1="repetition_miss20maf0.01",ARG2="filtered_vcf_miss20_maf0.01",ARG3="miss20.maf0.01" run_pixy_chr12_v02.sh
sbatch --export=ALL,ARG1="repetition_miss20maf0.01",ARG2="filtered_vcf_miss20_maf0.01",ARG3="miss20.maf0.01" run_pixy_chr17_v02.sh
sbatch --export=ALL,ARG1="repetition_miss20maf0.01",ARG2="filtered_vcf_miss20_maf0.01",ARG3="miss20.maf0.01" run_pixy_chr6_v02.sh
sbatch --export=ALL,ARG1="repetition_miss20maf0.01",ARG2="filtered_vcf_miss20_maf0.01",ARG3="miss20.maf0.01" run_pixy_chr23_v02.sh
~~~

Inside each folder containing the output files, I run this R script to concatenate the files per chromosome.

### Fst between Atlantic herring populations

run_all_herring_populations.sh
~~~bash
#!/bin/bash -l
 
## Load required modules
ml load bioinfo-tools pixy/1.2.5.beta1
## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Read arguments from the command line
DIR1=$ARG1
DIR2=$ARG2
FILTER=$ARG3

# Define directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}
PopDIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/pixy/"${DIR1}"/population_files"
input_vcf_dir="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/"${DIR2}

# What are the arguments?
echo ${WD}
echo ${PopDIR}
echo ${input_vcf_dir}

cd ${WD}

pixy --stats fst --populations ${PopDIR}/61_individuals.txt --vcf ${input_vcf_dir}/herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.${ChrName}.${ARG3}.vcf.gz --window_size 20000 --n_cores 4 --output_folder all_herring_populations --output_prefix all_herring_populations.${ChrName}.20kb.popgenpixy.out
~~~


### Merge tables

`wg_tables.R`
~~~R
#! /sw/apps/R/4.1.1/rackham/bin/Rscript
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

prefix <- args[1]
outname <- args[2]

files<-list.files(pattern=prefix)
tables<-lapply(files,read.table, header=T)
dxy<-rbindlist(tables,use.names=TRUE)
write.table(dxy,file=outname,col.names = T,row.names=F,quote=F,sep="\t")
~~~

This is how you run the script:

~~~bash
wg_tables.R 20kb.popgenpixy.out_pi.txt ${DIR2}.${DIR1}.wg.20kb.popgenpixy.out_pi.txt
wg_tables.R 20kb.popgenpixy.out_dxy.txt ${DIR2}.${DIR1}.wg.20kb.popgenpixy.out_dxy.txt
wg_tables.R 20kb.popgenpixy.out_fst.txt ${DIR2}.${DIR1}.wg.20kb.popgenpixy.out_fst.txt
~~~

### Plotting

