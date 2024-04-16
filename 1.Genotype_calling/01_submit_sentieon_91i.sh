#!/bin/bash
#SBATCH -A snic2020-5-94
#SBATCH -p core -n 10
#SBATCH -M rackham
#SBATCH -t 1-00:00:00
#SBATCH -J sentieon
#SBATCH -e sentieon_%A_%a.err
#SBATCH -o sentieon_%A_%a.out

ml python/2.7.15

TOPDIR=/crex/proj/snic2020-2-19/private/herring/users/mafalda/variants/91_indv
cd $TOPDIR

#V1.0 to match other samples
REFERENCE=/crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta
realpath /crex/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/clean_short_reads/mapped_untrim_reads/*bam > herring_bam.list
FILENAME=`cat $TOPDIR/herring_bam.list | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

SAMPBASE=$(basename $FILENAME)
SAMPLE=${SAMPBASE%.MD.RG.bam*}
SAMPLE_DIR=$(dirname $FILENAME)

echo $SAMPBASE
echo $SAMPLE
echo $SAMPLE_DIR

#should have ran this
sh /crex/proj/snic2020-2-19/private/herring/users/mafalda/software/Sentieon_pipeline_91indvs/herring_sentieon_91i.sh $FILENAME $REFERENCE $SAMPLE $SAMPBASE

