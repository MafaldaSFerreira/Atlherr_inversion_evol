#!/bin/bash
#SBATCH -A snic2020-5-94
#SBATCH -p node
#SBATCH -t 5-00:00:00
#SBATCH -J genotyper
#SBATCH -e genotyper_%A_%a.err            # File to which STDERR will be written
#SBATCH -o genotyper_%A_%a.out

ml python/2.7.15

TOPDIR=/crex/proj/snic2020-2-19/private/herring/users/mafalda/variants/91_indv
cd $TOPDIR

REFERENCE=/crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

#sentieon options
export SENTIEON_LICENSE=/domus/h1/mafaldaf/sentieon_files/Uppsala_University_Leif_Andersson_Lab_cluster3.lic
SENTIEON_INSTALL_DIR=/domus/h1/mafaldaf/sentieon_files/sentieon-genomics-201911
SENTIEON_TMPDIR=$SNIC_TMP

#make list of variants. But only do once
# for GVCF in */*.g.vcf.gz
# do
#   LST=gvcf_count_V2.txt
#   ls $GVCF >> gvcf_count_V2.txt
#   echo -ne " -v $GVCF" >> gvcf_list_V2.txt
# done

# wc -l gvcf_count_V2.txt

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REFERENCE -t 20 --algo GVCFtyper --emit_mode CONFIDENT herring_sentieon_91ind_190521.vcf.gz $(<gvcf_list_V2.txt)
