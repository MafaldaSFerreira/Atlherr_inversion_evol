# Genotype calling


### Mapping Celtic Sea and Baltic Sea individuals 

Mapping new data from 12 individuals from Celtic Sea and Baltic Sea (CS and BS individuals)

~~~bash
bwa/0.7.17
samtools/1.10

ml load bioinfo-tools
ml load bwa
ml load samtools
ml load gnuparallel

parallel -j 4 --colsep="\t" "bwa mem -M -t 4 -R ""{2}"" /proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta {1}_R1.fastq.gz {1}_R2.fastq.gz | samtools view -bS - > mapped_untrim_reads/{1}_pe.bam; samtools sort mapped_untrim_reads/{1}_pe.bam -o mapped_untrim_reads/{1}_pe_sorted.bam" :::: read_groups.txt
~~~

Mark duplicates

~~~bash
ml load picard/2.10.3
ml load gnuparallel

parallel -j 5 '$i={1}; java -jar $PICARD_ROOT/picard.jar MarkDuplicates INPUT=$i OUTPUT=${i/_pe_sorted.bam/.RG.MD.bam} METRICS_FILE=${i/_pe_sorted.bam/.metrics.file} CREATE_INDEX=true MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 TAG_DUPLICATE_SET_MEMBERS=false REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=DontTag REMOVE_DUPLICATES=false ASSUME_SORTED=false DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_MD5_FILE=false' ::: $(ls *_sorted.bam)
~~~

### Genotype calling using sentieon

We will use the bam files from [Han et al 2020](https://elifesciences.org/articles/61076). 

Put bam files in the same folder:
~~~bash
ln -s /proj/snic2020-2-19/private/herring/alignment/79_individuals/*.bam ./
ln -s /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/clean_short_reads/mapped_untrim_reads/*.bam .
~~~

Run Sentieon

First step. This first step calls HC Variant caller
~~~bash
sbatch --array=1-12 01_submit_sentieon_91i.sh
~~~

Then we need to run the joint call:

~~~bash
sbatch 02_genotyper_91.sh
~~~

And finally we run the hard filters using GATK:

~~~bash
sbatch 03_set_filter_91i_filter2.sh
~~~

We will run a final filter where we set the maximum (3 times the average depth for each sample) and minimum depth for each sample and remove indels

Extract depth for each sample
~~~bash
bcftools stats -s - herring_sentieon_91ind_190521.SV.VF.F2.setGT.inv.vcf.gz > herring_sentieon_91ind_190521.SV.VF.F2.setGT.inv.stats
~~~

Split vcf file for chromosomes:

~~~bash
ml load bioinfo-tools bcftools/1.12
ml load gnuparallel

parallel -j 20 'bcftools view -O z -o chromosomes/herring_sentieon_91ind_190521.SV.VF.F2.setGT.{1}.vcf.gz -r {1} herring_sentieon_91ind_190521.SV.VF.F2.setGT.vcf.gz' :::: chromosomes_26.txt
~~~

Run a loop to filter individuals
~~~bash
for i in $(ls *.vcf.gz); do bcftools filter -S . -e "FMT/DP[0] > 30.9" $i | bcftools filter -S . -e "FMT/DP[1] > 36.3" | bcftools filter -S . -e "FMT/DP[2] > 36" | bcftools filter -S . -e "FMT/DP[3] > 36.6" | bcftools filter -S . -e "FMT/DP[4] > 34.2" | bcftools filter -S . -e "FMT/DP[5] > 34.5" | bcftools filter -S . -e "FMT/DP[6] > 30.3" | bcftools filter -S . -e "FMT/DP[7] > 33" | bcftools filter -S . -e "FMT/DP[8] > 39" | bcftools filter -S . -e "FMT/DP[9] > 34.2" | bcftools filter -S . -e "FMT/DP[10] > 38.4" | bcftools filter -S . -e "FMT/DP[11] > 32.4" | bcftools filter -S . -e "FMT/DP[12] > 30.9" | bcftools filter -S . -e "FMT/DP[13] > 32.1" | bcftools filter -S . -e "FMT/DP[14] > 30.6" | bcftools filter -S . -e "FMT/DP[15] > 39.9" | bcftools filter -S . -e "FMT/DP[16] > 33.9" | bcftools filter -S . -e "FMT/DP[17] > 33" | bcftools filter -S . -e "FMT/DP[18] > 30.9" | bcftools filter -S . -e "FMT/DP[19] > 35.1" | bcftools filter -S . -e "FMT/DP[20] > 34.5" | bcftools filter -S . -e "FMT/DP[21] > 32.1" | bcftools filter -S . -e "FMT/DP[22] > 91.5" | bcftools filter -S . -e "FMT/DP[23] > 84.6" | bcftools filter -S . -e "FMT/DP[24] > 87.3" | bcftools filter -S . -e "FMT/DP[25] > 94.8" | bcftools filter -S . -e "FMT/DP[26] > 75.6" | bcftools filter -S . -e "FMT/DP[27] > 84.3" | bcftools filter -S . -e "FMT/DP[28] > 79.5" | bcftools filter -S . -e "FMT/DP[29] > 86.7" | bcftools filter -S . -e "FMT/DP[30] > 90.3" | bcftools filter -S . -e "FMT/DP[31] > 130.5" | bcftools filter -S . -e "FMT/DP[32] > 129.6" | bcftools filter -S . -e "FMT/DP[33] > 88.8" | bcftools filter -S . -e "FMT/DP[34] > 126.6" | bcftools filter -S . -e "FMT/DP[35] > 78.6" | bcftools filter -S . -e "FMT/DP[36] > 144" | bcftools filter -S . -e "FMT/DP[37] > 90.9" | bcftools filter -S . -e "FMT/DP[38] > 96" | bcftools filter -S . -e "FMT/DP[39] > 101.1" | bcftools filter -S . -e "FMT/DP[40] > 69.6" | bcftools filter -S . -e "FMT/DP[41] > 79.5" | bcftools filter -S . -e "FMT/DP[42] > 80.1" | bcftools filter -S . -e "FMT/DP[43] > 138.6" | bcftools filter -S . -e "FMT/DP[44] > 75.6" | bcftools filter -S . -e "FMT/DP[45] > 115.5" | bcftools filter -S . -e "FMT/DP[46] > 72.9" | bcftools filter -S . -e "FMT/DP[47] > 99.9" | bcftools filter -S . -e "FMT/DP[48] > 88.2" | bcftools filter -S . -e "FMT/DP[49] > 105.6" | bcftools filter -S . -e "FMT/DP[50] > 99" | bcftools filter -S . -e "FMT/DP[51] > 106.2" | bcftools filter -S . -e "FMT/DP[52] > 97.8" | bcftools filter -S . -e "FMT/DP[53] > 87.3" | bcftools filter -S . -e "FMT/DP[54] > 112.5" | bcftools filter -S . -e "FMT/DP[55] > 99.6" | bcftools filter -S . -e "FMT/DP[56] > 129" | bcftools filter -S . -e "FMT/DP[57] > 112.8" | bcftools filter -S . -e "FMT/DP[58] > 117.9" | bcftools filter -S . -e "FMT/DP[59] > 111.9" | bcftools filter -S . -e "FMT/DP[60] > 102.6" | bcftools filter -S . -e "FMT/DP[61] > 113.7" | bcftools filter -S . -e "FMT/DP[62] > 122.1" | bcftools filter -S . -e "FMT/DP[63] > 100.5" | bcftools filter -S . -e "FMT/DP[64] > 107.7" | bcftools filter -S . -e "FMT/DP[65] > 111.6" | bcftools filter -S . -e "FMT/DP[66] > 93.6" | bcftools filter -S . -e "FMT/DP[67] > 103.8" | bcftools filter -S . -e "FMT/DP[68] > 109.2" | bcftools filter -S . -e "FMT/DP[69] > 106.2" | bcftools filter -S . -e "FMT/DP[70] > 48.6" | bcftools filter -S . -e "FMT/DP[71] > 49.8" | bcftools filter -S . -e "FMT/DP[72] > 49.5" | bcftools filter -S . -e "FMT/DP[73] > 30.6" | bcftools filter -S . -e "FMT/DP[74] > 38.4" | bcftools filter -S . -e "FMT/DP[75] > 30.3" | bcftools filter -S . -e "FMT/DP[76] > 118.5" | bcftools filter -S . -e "FMT/DP[77] > 57.6" | bcftools filter -S . -e "FMT/DP[78] > 99.6" | bcftools filter -S . -e "FMT/DP[79] > 102.9" | bcftools filter -S . -e "FMT/DP[80] > 119.1" | bcftools filter -S . -e "FMT/DP[81] > 123" | bcftools filter -S . -e "FMT/DP[82] > 110.4" | bcftools filter -S . -e "FMT/DP[83] > 34.8" | bcftools filter -S . -e "FMT/DP[84] > 120.6" | bcftools filter -S . -e "FMT/DP[85] > 90.6" | bcftools filter -S . -e "FMT/DP[86] > 102.6" | bcftools filter -S . -e "FMT/DP[87] > 97.5" | bcftools filter -S . -e "FMT/DP[88] > 33" | bcftools filter -S . -e "FMT/DP[89] > 36.9" | bcftools filter -S . -e "FMT/DP[90] > 51" | bcftools filter -S . -e "FMT/DP < 3" | bcftools view -f PASS -e 'ALT="*" | TYPE~"indel" | ref="N"' -O z -o ${i/herring_sentieon_91ind_190521.SV.VF.F2.setGT./herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.}; done
~~~

After, perform a filtering based on missing data and minor allele frequency for downstream analyses. We include variant and invariant sites:

#### maf 0.1
~~~bash
#!/bin/bash -l

ml load bioinfo-tools vcftools/0.1.16 bcftools/1.12

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs"
outdir_WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/filtered_vcf_miss20maf0.1"

ChrName=chr${SLURM_ARRAY_TASK_ID}

cd ${WD}

vcf="herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv."${ChrName}".vcf.gz"
vcf_inv="herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv."${ChrName}".inv.vcf.gz"
vcf_var="herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv."${ChrName}".var.vcf.gz"
vcf_output="herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv."${ChrName}".miss20.maf0.05.vcf.gz"

echo ${WD}
echo ${outdir_WD}
echo ${vcf} 
echo ${vcf_inv}
echo ${vcf_var}
echo ${vcf_output}

vcftools --gzvcf ${vcf} --max-maf 0 --max-missing 0.8 --recode --stdout | bgzip -c > ${outdir_WD}/${vcf_inv}
vcftools --gzvcf ${vcf} --mac 1 --max-missing 0.8 --maf 0.1 --recode --stdout | bgzip -c > ${outdir_WD}/${vcf_var}

tabix ${outdir_WD}/${vcf_inv}
tabix ${outdir_WD}/${vcf_var}

# combine the two VCFs using bcftools concat
bcftools concat -O z -o ${outdir_WD}/${vcf_output} --allow-overlaps ${outdir_WD}/${vcf_inv} ${outdir_WD}/${vcf_var} 

tabix ${outdir_WD}/${vcf_output}

#rm ${outdir_WD}/${vcf_inv}
#rm ${outdir_WD}/${vcf_var}
~~~

#### maf 0.01
~~~bash
#!/bin/bash -l

ml load bioinfo-tools vcftools/0.1.16 bcftools/1.12

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs"
outdir_WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/GenomicsGenerals/F2_vcfs/filtered_vcf_miss20_maf0.01"

ChrName=chr${SLURM_ARRAY_TASK_ID}

cd ${WD}

vcf="herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv."${ChrName}".vcf.gz"
vcf_inv="herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv."${ChrName}".inv.vcf.gz"
vcf_var="herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv."${ChrName}".var.vcf.gz"
vcf_output="herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv."${ChrName}".miss20.maf0.01.vcf.gz"

echo ${WD}
echo ${outdir_WD}
echo ${vcf} 
echo ${vcf_inv}
echo ${vcf_var}
echo ${vcf_output}

vcftools --gzvcf ${vcf} --max-maf 0 --max-missing 0.8 --recode --stdout | bgzip -c > ${outdir_WD}/${vcf_inv}

vcftools --gzvcf ${vcf} --mac 1 --max-missing 0.8 --maf 0.01 --recode --stdout | bgzip -c > ${outdir_WD}/${vcf_var}

tabix ${outdir_WD}/${vcf_inv}
tabix ${outdir_WD}/${vcf_var}

# combine the two VCFs using bcftools concat
bcftools concat --allow-overlaps ${outdir_WD}/${vcf_inv} ${outdir_WD}/${vcf_var} -O z -o ${outdir_WD}/${vcf_output}

tabix ${outdir_WD}/${vcf_output}

rm ${outdir_WD}/${vcf_inv}
rm ${outdir_WD}/${vcf_var}

bcftools stats ${outdir_WD}/${vcf_output} > ${outdir_WD}/${vcf_output/.vcf.gz/.stats.txt}
~~~


### Call consensus fasta files

`create_consensus.sh` is in this folder.

the script uses `do_bed.awk` which can be found in in my [github](https://github.com/MafaldaSFerreira/wtjr_winter_camouflage_evolution/blob/master/variant_call_and_consensus_fasta/do_bed.awk) .

~~~bash
ml load bioinfo-tools bcftools/1.12
ml load gnuparallel

parallel -j 20 'bcftools index herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.{1}.vcf.gz' :::: chromosomes_26.txt

parallel -j 20 'bcftools view -s {1} --exclude-uncalled -O z -o {2}/{1}.{2}.vcf.gz herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.{2}.vcf.gz' :::: individuals.txt :::: chromosomes_26.txt

parallel -j 1 "/proj/snic2020-2-19/private/herring/users/mafalda/software/create_consensus.sh {1}" :::: chromosomes_26.txt
~~~

