#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
set -x
data_dir=/crex/proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/clean_short_reads/mapped_untrim_reads
# MF: I dont know what this dir is, but assuming is the output, I'll give 
# a potential output.
TOPDIR=/crex/proj/snic2020-2-19/private/herring/users/mafalda/variants/91_indv

#fastq_1=$1
#fastq_2=$2 #If using Illumina paired data

# Update with the location of the reference data files
fasta=$2

# Set SENTIEON_LICENSE if it is not set in the environment
export SENTIEON_LICENSE=/domus/h1/mafaldaf/sentieon_files/Uppsala_University_Leif_Andersson_Lab_cluster3.lic

# Update with the location of the Sentieon software package
SENTIEON_INSTALL_DIR=/domus/h1/mafaldaf/sentieon_files/sentieon-genomics-201911

# Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=$SNIC_TMP

# It is important to assign meaningful names in actual cases.
# It is particularly important to assign different read group names.
sample=$3
group=$4
platform="ILLUMINA"

# Other settings
nt=10 #number of threads to use in computation

# ******************************************
# 0. Setup
# ******************************************
workdir=$TOPDIR/$sample #ede edited this to be sample
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

#Sentieon proprietary compression
bam_option="--bam_compression 1"

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000

# speed up memory allocation malloc in bwa
export LD_PRELOAD=$SENTIEON_INSTALL_DIR/lib/libjemalloc.so
export MALLOC_CONF=lg_dirty_mult:-1
#
# ( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_1 $fastq_2 || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort $bam_option -r $fasta -o ${sample}.sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2. Metrics
# ******************************************
# $SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i ${sample}.sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
# $SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
# $SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
# $SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
# $SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads
# To mark duplicate reads only without removing them, remove "--rmdup" in the second command
# ******************************************
# $SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i ${sample}.sorted.bam --algo LocusCollector --fun score_info score.txt
# $SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i ${sample}.sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt $bam_option ${sample}.MD.RG.bam

# ******************************************
# 6. HC Variant caller
# Note: Sentieon default setting matches versions before GATK 3.7.
# Starting GATK v3.7, the default settings have been updated multiple times.
# Below shows commands to match GATK v3.7 - 4.1
# Please change according to your desired behavior.
# ******************************************

# Matching GATK 3.7, 3.8, 4.0
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i ${sample}.MD.RG.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=10 output-hc.vcf.gz

# Matching GATK 4.1
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $data_dir/${sample}.MD.RG.bam --algo Haplotyper --genotype_model multinomial --emit_mode gvcf --emit_conf 30 --call_conf 30 ${sample}.g.vcf.gz
