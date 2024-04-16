#!/bin/bash
#SBATCH -A snic2020-5-94
#SBATCH -p core -n 4 #was 4
#SBATCH -t 5-00:00:00
#SBATCH -J set_filters
#SBATCH -e set_filters_%J_%A_%a.err
#SBATCH -o set_filters_%J_%A_%a.out

module load bioinfo-tools samtools vcftools bcftools picard/2.20.4
module load GATK/4.1.4.1

TOPDIR=/crex/proj/snic2020-2-19/private/herring/users/mafalda/variants/91_indv
cd $TOPDIR
REF=/crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

gatk --java-options "-Xmx29g" VariantFiltration  \
    -R $REF \
    -V herring_sentieon_91ind_190521.SV.vcf.gz \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0))" \
    --filter-name "RPRS8" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('QD') && QD < 2.0))" \
    --filter-name "QD2" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('FS') && FS > 60.0))" \
    --filter-name "FS60" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('SOR') && SOR > 3.0))" \
    --filter-name "SOR3" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('MQ') && MQ < 40.0))" \
    --filter-name "MQ40" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
    --filter-name "MQ12.5" \
    -G-filter "vc.isSNP() && DP < 3" \
    -G-filter-name "gtDP1" \
    -G-filter "vc.isSNP() && GQ < 20" \
    -G-filter-name "gtGQ10" \
    -O herring_sentieon_91ind_190521.SV.VF.F2.vcf.gz

#set filtered GT to no call
gatk --java-options "-Xmx29g" SelectVariants  \
    -R $REF \
    -V  herring_sentieon_91ind_190521.SV.VF.F2.vcf.gz \
    --set-filtered-gt-to-nocall \
    -O herring_sentieon_91ind_190521.SV.VF.F2.setGT.vcf.gz
