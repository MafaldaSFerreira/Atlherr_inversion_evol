#for chromosome as $1
chromosome=$1
#parallel 'echo "creating wwe consensus input for {}"' ::: $chromosome

cd "$chromosome"
echo $chromosome
echo "I'm inside" 

ml load bioinfo-tools bcftools/1.12 BEDTools/2.29.2 samtools/1.12
ml load gnuparallel

# create the bed file with all positions with a call:
parallel -j 20 'i={1}; bcftools view -H $i | /proj/snic2020-2-19/private/herring/users/mafalda/software/wtjr_camouflage/variant_call_and_consensus_fasta/do_bed.awk > ${i/.vcf.gz/.filter.bed}' ::: $(ls *.vcf.gz)

# grep the length of chromosome:
grep -w "$chromosome" /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/Ch_v2.0.2.chromsizes > "$chromosome".chromsize

# use bedtools complement to create a bed file of positions that had no call:
#echo "creating complement bed" 
parallel -j 20 --colsep="\t" 'i={1} ; bedtools complement -i $i -g {2}.chromsize > ${i/.bed/.nocall.bed}' ::: $(ls *.filter.bed) ::: $chromosome

parallel -j 20 'bcftools index {1}' ::: $(ls *.vcf.gz)

# call the consensus fasta:
# individuals.txt is a file with individual names
parallel -j 20 --colsep="\t" 'samtools faidx /crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta {2} | bcftools consensus -I {1}.{2}.vcf.gz --sample {1} > {1}.{2}.fa' :::: ../individuals.txt ::: $chromosome

parallel -j 20 --colsep="\t" 'i={1}.{2}.fa ; bedtools maskfasta -fi $i -bed ${i/.fa/.filter.nocall.bed} -fo ${i/.fa/.consensus.fa}' :::: ../individuals.txt ::: $chromosome

# remove bed files because they are quite bid.
parallel -j 20 --colsep="\t" 'rm {1}.{2}.fa' :::: ../individuals.txt ::: $chromosome
rm *.bed
rm *.gz.csi
# samtools works with bgzipped files:
for i in $(ls *.consensus.fa); do bgzip $i; done