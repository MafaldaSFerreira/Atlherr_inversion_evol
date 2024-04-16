

import sys
import pandas as pd

# colnames for blast result format 6 table
colnames=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
 "send", "evalue", "bitscore"]

# read blast results
blast_result_start_hap1=pd.read_csv(sys.argv[1], names=colnames, header=None, sep="\t")
blast_result_end_hap1=pd.read_csv(sys.argv[2], names=colnames, header=None, sep="\t")

blast_result_start_hap2=pd.read_csv(sys.argv[3], names=colnames, header=None, sep="\t")
blast_result_end_hap2=pd.read_csv(sys.argv[4], names=colnames, header=None, sep="\t")

filename=sys.argv[1]
names=filename.split(".")

# HAP 1 ####
## START ####
## sort by send, evalue and bitscore. 
## in certain cases, if there are two equal send, this means we will be chosing the best hit, 
## even if the second best includes a coordinate that would result in a longer sequence.
## I think it's probably best to be conservative.

blast_sorted_start_hap1=blast_result_start_hap1.sort_values(by=['sstart', 'evalue', 'bitscore'], ascending=[True, True, False])

## select the coordinate of the first entry of the table:
coordinate_start_hap1=blast_sorted_start_hap1.loc[:,'qstart'].iloc[0]


## END
## sort by send, evalue and bitscore. 
## in certain cases, if there are two equal send, this means we will be chosing the best hit, 
## even if the second best includes a coordinate that would result in a longer sequence.
## I think it's probably best to be conservative.

blast_sorted_end_hap1=blast_result_end_hap1.sort_values(by=['send', 'evalue', 'bitscore'], ascending=[False, True, False])

## select the coordinate of the first entry of the table:
coordinate_end_hap1=blast_sorted_end_hap1.loc[:,'qend'].iloc[0]

## in the case that we have a coordinate end that comes before the coordinate start, we just take the best hit using the bitscore, and the the evalue, as is originally outut by blast
if coordinate_end_hap1 < coordinate_start_hap1:
    blast_sorted_end_hap1=blast_result_end_hap1.sort_values(by=['bitscore', 'evalue'], ascending=[False, True])
    
    coordinate_end_hap1=blast_sorted_end_hap1.loc[:,'qend'].iloc[0]


# HAP 2 ####
## START ####

blast_sorted_start_hap2=blast_result_start_hap2.sort_values(by=['sstart', 'evalue', 'bitscore'], ascending=[True, True, False])

## select the coordinate of the first entry of the table:
coordinate_start_hap2=blast_sorted_start_hap2.loc[:,'qstart'].iloc[0]


## END
## sort by send, evalue and bitscore. 
## in certain cases, if there are two equal send, this means we will be chosing the best hit, 
## even if the second best includes a coordinate that would result in a longer sequence.
## I think it's probably best to be conservative.

blast_sorted_end_hap2=blast_result_end_hap2.sort_values(by=['send', 'evalue', 'bitscore'], ascending=[False, True, False])

## select the coordinate of the first entry of the table:
coordinate_end_hap2=blast_sorted_end_hap2.loc[:,'qend'].iloc[0]

## in the case that we have a coordinate end that comes before the coordinate start, we just take the best hit using the bitscore, and the the evalue, as is originally outut by blast
if coordinate_end_hap2 < coordinate_start_hap2:
    blast_sorted_end_hap2=blast_result_end_hap2.sort_values(by=['bitscore', 'evalue'], ascending=[False, True])
    
    coordinate_end_hap2=blast_sorted_end_hap2.loc[:,'qend'].iloc[0]

# Calculate the length of each haplotype
length_hap1=coordinate_end_hap1-coordinate_start_hap1
length_hap2=coordinate_end_hap2-coordinate_start_hap2

if length_hap1 > length_hap2:

    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(names[0], "hap1", coordinate_start_hap1, coordinate_end_hap1, length_hap1,
    "hap2", coordinate_start_hap2, coordinate_end_hap2, length_hap2))

else:

    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(names[0], "hap2", coordinate_start_hap2, coordinate_end_hap2, length_hap2,
    "hap1", coordinate_start_hap1, coordinate_end_hap1, length_hap1))

