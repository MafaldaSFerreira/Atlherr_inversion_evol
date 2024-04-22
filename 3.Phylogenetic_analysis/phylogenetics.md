# Phylogenetics 

Concatenate whole genome chromosome sequences for each individual. These sequences where generated in `2ii.ancestral_state_from_sprat.md` and `genotype_call_pipeline.md`.

The headers of the fasta files contain the name for each individual.

~~~bash
## Load required modules
ml load bioinfo-tools bcftools/1.12 BEDTools/2.29.2 samtools/1.12

## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus"
output_DIR=${WD}"/full_alignments_27022023"
input_DIR=${WD}"/"${ChrName}

# Generate alignments without sprat
cat ${input_DIR}/*.named.fa.gz > ${output_DIR}/${ChrName}".91indvs.algn.fa.gz"

# Generate alignments with sprat

cat ${input_DIR}/*.named.fa.gz ${WD}/Sprat/EuSprat.${ChrName}.named.fa.gz > ${output_DIR}/${ChrName}".91indvs.outgroup.algn.fa.gz"
~~~

Unzip the files to run AMAS

~~~bash
for i in $(ls *fa.gz); do bgzip -d $i; done
~~~

### Whole genome alignments

Trim all missing data from the whole-genome alignments, with and without the outgroup (sprat).

~~~bash
#!/bin/bash -l
 
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/full_alignments_27022023"
cd ${WD}

/proj/snic2020-2-19/private/herring/users/mafalda/software/AMAS/amas/AMAS.py concat -i *.no_outgroup.fa -f fasta -d dna -p whole_genome_alignment.91indvs.no_outgroup.nomiss.txt -t whole_genome_alignment.91indvs.no_outgroup.nomiss.fa -u fasta -y unspecified -c 18

/proj/snic2020-2-19/private/herring/users/mafalda/software/AMAS/amas/AMAS.py concat -i *.outgroup.fa -f fasta -d dna -p whole_genome_alignment.91indvs.outgroup.nomiss.txt -t whole_genome_alignment.91indvs.outgroup.nomiss.fa -u fasta -y unspecified -c 18
~~~

### Inversion alignments

run_create_inversion_alignments.sh
~~~bash
#!/bin/bash -l
 
ml load bioinfo-tools bcftools/1.12 samtools/1.12

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus"
OUTGROUP_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/Sprat"
OUTPUT_DIR_INDV=${WD}"/inversion_alignments/individual_alignments"
OUTPUT_DIR_FULL=${WD}"/inversion_alignments/full_alignments"

chromosomes=("chr6" "chr12" "chr17" "chr23")
coordinates=(6 12 17 23)
start=(22282765 17826318 25805445 16226443)
end=(24868581 25603093 27568510 17604273)

chrName=${chromosomes[$SLURM_ARRAY_TASK_ID]} 
c=${coordinates[$SLURM_ARRAY_TASK_ID]} 
s=${start[$SLURM_ARRAY_TASK_ID]}
e=${end[$SLURM_ARRAY_TASK_ID]}

# Chromosome specific directories
Chr_DIR=${WD}/${chrName}
OUTPUT_DIR_INDV_CHR=${OUTPUT_DIR_INDV}/${chrName}

# Retrieve the alignemnts of the inversions
cd ${Chr_DIR}

for i in $(ls *.modified.fa.gz); do o=${i/.modified.fa.gz/.tmp.fa}; samtools faidx $i ${c}:${s}-${e} > ${OUTPUT_DIR_INDV_CHR}/$o; done

# Modify the headers 
cd ${OUTPUT_DIR_INDV_CHR}

for i in $(ls *.tmp.fa); do header=${i%.${chrName}.tmp*}; o=${i/.tmp.fa/.TMP.fa}; sed "s/>.*/>$header/" $i > $o; done

#Concatenate the alignments
cat *TMP.fa > ${OUTPUT_DIR_FULL}/${chrName}.inversion.alignment.no_outgroup.fa

# Add the outgroup

cd ${OUTGROUP_DIR}

samtools faidx EuSprat.${chrName}.modified.fa.gz ${c}:${s}-${e} > ${OUTPUT_DIR_INDV_CHR}/EuSprat.${chrName}.tmp.fa

cd ${OUTPUT_DIR_INDV_CHR}

sed "s/>.*/>EuSprat/" EuSprat.${chrName}.tmp.fa > EuSprat.${chrName}.TMP.fa

cat EuSprat.${chrName}.TMP.fa ${OUTPUT_DIR_FULL}/${chrName}.inversion.alignment.no_outgroup.fa > ${OUTPUT_DIR_FULL}/${chrName}.inversion.alignment.outgroup.fa
~~~

Filter inversion alignments. We allow 50% missing data.

run_filter_inv_algns.sh
~~~bash
#!/bin/bash -l
 
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/inversion_alignments/full_alignments"
OUTDIR=${WD}"/filtered_alignments"

cd ${WD}

files=($(ls *.fa))

input_alignment=${files[$SLURM_ARRAY_TASK_ID]}

/proj/snic2020-2-19/private/herring/users/mafalda/software/AMAS/amas/AMAS.py trim -i ${input_alignment} -u fasta -o ${OUTDIR}/${input_alignment/.fa/.miss50.fa} -t 0.5 -f fasta -d dna

~~~

### Generate Maximum Likelihood trees with iqtree

Run inversions without outgroup:

~~~bash
#!/bin/bash -l

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/inversion_alignments/ML_trees_Pacific_outgroup"
INPUT_DIR_1="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/inversion_alignments/full_alignments/filtered_alignments"
INPUT_DIR_2="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/inversion_alignments/full_alignments"

cd ${WD}

files=($(ls ${INPUT_DIR_1}/*.no_outgroup.miss50.fa*) $(ls ${INPUT_DIR_2}/*.no_outgroup.*))

input_alignment=${files[$SLURM_ARRAY_TASK_ID]}

iqtree -s ${input_alignment} -m MFP -nt AUTO -bb 100 -redo -o Pacific5_Vancouver_Pacific
~~~

Run inversions with outgroup:

~~~bash
#!/bin/bash -l
 
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/inversion_alignments/ML_trees_Sprat_outgroup"
INPUT_DIR_1="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/inversion_alignments/full_alignments/filtered_alignments"
INPUT_DIR_2="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/inversion_alignments/full_alignments"

cd ${WD}

files=($(ls ${INPUT_DIR_1}/*.outgroup.miss50.fa*) $(ls ${INPUT_DIR_2}/*.outgroup.*))

input_alignment=${files[$SLURM_ARRAY_TASK_ID]}

iqtree -s ${input_alignment} -m MFP -nt AUTO -bb 100 -redo -o EuSprat
~~~

Run whole-genome alignments without outgroup (we specificy one Pacific herring individual as outgroup)

~~~bash
#!/bin/bash -l

ml load bioinfo-tools iqtree/2.0-rc2-omp-mpi

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/full_alignments_27022023/whole_genome_alignments"

cd ${WD}

iqtree -s whole_genome_alignment.91indvs.no_outgroup.nomiss.fa -m MFP -nt AUTO -bb 1000 -redo -o Pacific5_Vancouver_Pacific -T 18
~~~

Run whole-genome alignments with outgroup:

~~~bash

#!/bin/bash -l
ml load bioinfo-tools iqtree/2.0-rc2-omp-mpi

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/full_alignments_27022023/whole_genome_alignments"

cd ${WD}

iqtree -s whole_genome_alignment.91indvs.outgroup.nomiss.fa -m MFP -nt AUTO -bb 1000 -redo -o EuSprat -T 18