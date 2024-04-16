#!/bin/env python3

# this script will take a CDS alignment with several individuals, generate a consensus for two specified
# populations and calculate dnds using an outgroup. here the outgroup is "EuSprat"
#
# Description of algorithm:
#
# 1. create a consensus sequence from a group of individuals that are defined as 
# a population. The consensus converts any position that contains IUPAC codes (hets or N), or that 
# contains alternative homozygote alleles (e.g. AAAATTTT) into an N. Simply, any position that is not fixed
# in the population will be discarded. This is quite conservative since even if one individual has 100% missing
# data, the entire gene will be discarded.
#
# 2. convert stop codons in the consensus into gaps "---"
# 
# 3. convert codons containing IUPAC codes in the consensus into gaps "---"
#
# 4. filter out any gaps from the alignment, preserving codons
#
# 5. calculate dNdS with the ML method (M0 in codeml) for alignments longer than 100 codons
#
#
# first arg: alignment file in fasta format
# second arg: individuals in population 1
# third arg: individuals in population 2
# fourth arg: output file
# Mafalda S. Ferreira, February 2022, Uppsala, SE

import sys
import Bio
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.codonalign import CodonAlignment, CodonSeq
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord


######## FUNCTIONS ########################################
# Function to convert stop codons into ambiguity codes
# took this code from here https://www.biostars.org/p/296261/
def replace_stop_codons(record, codon_stop_array = ["TAG", "TGA", "TAA"]):
    tempRecordSeq = list(record.seq)
    for index in range(0, len(record.seq), 3):
            codon = record.seq[index:index+3]
            if codon in codon_stop_array:
                tempRecordSeq[index:index+3] = '-','-','-'
    record.seq = Seq("".join(tempRecordSeq))
    return record

def replace_IUPAC(record, IUPAC_array = ["R","Y","S","W","K","M","B","D","H","V","N"]):
    tempRecordSeq = list(record.seq)
    for index in range(0, len(record.seq), 3):
            codon = record.seq[index:index+3]
            for n in range(0, len(codon)):
                base = codon[n]
                if base in IUPAC_array:
                    tempRecordSeq[index:index+3] = '-','-','-'
    record.seq = Seq("".join(tempRecordSeq))
    return record
  
######## INPUTS ########################################
# Read input file
inputfile = sys.argv[1]
pop1 = open(sys.argv[2], "r")
pop2 = open(sys.argv[3], "r")

pop1_array = pop1.read().split()
pop2_array = pop2.read().split()

outputfile = open(sys.argv[4],"w")

# Read input file as a SeqRecord
#alignment=SeqIO.parse(inputfile, "fasta")
alignment = AlignIO.read(inputfile, "fasta")

######## CREATE CONSENSUS ########################################
# Create consensus:
mapping = {}
index=0

for record in alignment:
    mapping[record.id] = index
    index += 1

#print(mapping)

## Pop1 consensus:
pop1_alignment = Bio.Align.MultipleSeqAlignment([])

for k in mapping.keys():
    if k in pop1_array:
        v=mapping[k]
        pop1_alignment.append(alignment[v])

#print(pop1_alignment)

summary_align = AlignInfo.SummaryInfo(pop1_alignment)
pop1_consensus = summary_align.dumb_consensus(threshold=1, ambiguous='N')

## Pop2 consensus:
pop2_alignment = Bio.Align.MultipleSeqAlignment([])

for k in mapping.keys():
    if k in pop2_array:
        v=mapping[k]
        pop2_alignment.append(alignment[v])

#print(pop2_alignment)

summary_align = AlignInfo.SummaryInfo(pop2_alignment)
pop2_consensus = summary_align.dumb_consensus(threshold=1, ambiguous='N')

## Outgroup consensus:
outgroup=["EuSprat"]

outgroup_alignment = Bio.Align.MultipleSeqAlignment([])

for k in mapping.keys():
    if k in outgroup:
        v=mapping[k]
        outgroup_alignment.append(alignment[v])

#print(outgroup_alignment)

# Create new alignment
#print(type(str(pop1_consensus)))
#print(type(str(pop2_consensus)))

pop1_record=SeqRecord(Seq(str(pop1_consensus)), id="pop1_consensus", name="pop1_consensus")
pop2_record=SeqRecord(Seq(str(pop2_consensus)), id="pop2_consensus", name="pop2_consensus")
outgroup_record=SeqRecord(Seq(str(outgroup_alignment[0].seq)), id="EuSprat", name="EuSprat")

dnds_align=[]

dnds_align.append(pop1_record)
dnds_align.append(pop2_record)
dnds_align.append(outgroup_record)

######## ADDITIONAL FILTERING ########################################
new_alignment = Bio.Align.MultipleSeqAlignment([])

# Remove stop_codons
# Another doubt I have is keeping or not IUPAC codes and missing data N. I know most people will 
# convert these to ambiguity codes and cleandata in paml, effectively ignoring these sites. 
# I think this is fair. 

for record in dnds_align:
    record_no_stops=replace_stop_codons(record)
    record_no_stops_no_IUPAC=replace_IUPAC(record_no_stops)
    #print(record_no_stops_no_IUPAC)
    new_alignment.append(record_no_stops_no_IUPAC)

gene=inputfile.split(".")[0]
output_consensus_name = "{}.consensus.fa".format(gene)
AlignIO.write(new_alignment, output_consensus_name, "fasta")

######## Check missing data #######################################

for record in new_alignment:
    # Take the current sequence
    sequence_name=str(record.id)
    sequence=str(record.seq)
    # Calculate the percentage of missing data for each sequence
    missing=float(sequence.count("-"))/float(len(sequence))*100
    print("Individual: {} Length: {} Missing: {}".format(sequence_name, len(sequence), missing))
    if(missing >= 50):
        print("There's a lot of missing data in sequence {}".format(sequence_name))

######## Filter missing codons #######################################

# Find codons with missing data and remove them from the alignment
# This code preserves codons in the alignment
MISS_array = ['-']
coord_to_remv = []
for record in new_alignment:
    for index in range(0, len(record.seq), 3):
            codon = record.seq[index:index+3]
            # set is a faster way to compare two lists and check if there are common elements 
            # https://stackoverflow.com/a/45099098
            if(set(list(codon)).intersection(set(MISS_array))):
                # print(codon)
                # store the coordinates that have gaps.
                coord_to_remv.append([index,index+3])

# This will remove duplicated lists within the lists
# https://stackoverflow.com/a/12198497
coord_to_remv_set = sorted(set(map(tuple,coord_to_remv)))  #need to convert the inner lists to tuples so they are hashable

# Solution from: https://www.biostars.org/p/419108/#419123
# Zipping each one will create two long lists with coordinates
# The end of gaps are the start of regions we want to keep
# And the start of gaps are the end of regions we want to keep
ends, starts = zip(*coord_to_remv_set)

# Create a new alignment with the beginning of the alignment
filtered_algn = new_alignment[:,:ends[0]]

for c in range(0,len(starts)-1):
    # append the new alignments:
    # for each start, the end of the region will be a member of the next tupple, so we add +1
    filtered_algn += new_alignment[:,starts[c]:ends[c+1]]

# Write the filtered alignment
output_consensus_filtered_name = "{}.consensus.filtered.fa".format(gene)
AlignIO.write(filtered_algn, output_consensus_filtered_name, "fasta")

# We should not run the analisis if the alignment is too short after filtering
if filtered_algn.get_alignment_length() <= 300:
    print("Alingment is less than 100 codons. Skipping")
    print("Alignment_length\t{} codons".format(filtered_algn.get_alignment_length()/3))
    
    pop1.close()
    pop2.close()
    outputfile.close()
else:

    ######## DnDs ratios #######################################
    # Create empty codon alignment object
    codon_aln = CodonAlignment()

    codon_aln_no_stops=codon_aln.from_msa(filtered_algn)

    dN_pop1, dS_pop1 = cal_dn_ds(codon_aln_no_stops[0], codon_aln_no_stops[2], method="ML")

    dNdS_pop1 = dN_pop1/dS_pop1

    dN_pop2, dS_pop2 = cal_dn_ds(codon_aln_no_stops[1], codon_aln_no_stops[2], method="ML")
    dNdS_pop2 = dN_pop2/dS_pop2

    output_header = "transcript\tdN_pop1\tdS_pop1\tdNdS_pop1\tdN_pop2\tdS_pop2\tdNdS_pop2\n"
    #output_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene, round(dN_pop1, ndigits=9), round(dS_pop1, ndigits=9), round(dNdS_pop1, ndigits=9), round(dN_pop2, ndigits=9), round(dS_pop2, ndigits=9), round(dNdS_pop2, ndigits=9))
    output_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene, dN_pop1, dS_pop1, dNdS_pop1, dN_pop2, dS_pop2, dNdS_pop2)

    outputfile.write(output_header)
    outputfile.write(output_line)

    ######## Close files #######################################

    pop1.close()
    pop2.close()
    outputfile.close()(mutation)