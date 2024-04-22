# A Python to conver an alignment from mafft into a vcf file
# in the coordinates of the Ch_v2.0.2 reference genome.
# You need to provide a file containing the coordinates of the 
# gene in case.
# Mafalda Ferreira

def split(word):
    return list(word)

import sys
from Bio import SeqIO, AlignIO
#import re
import pandas as pd

input_alig=open(sys.argv[1],"r")
input_info=open(sys.argv[2],"r")
output=open(sys.argv[3],'w')

align = AlignIO.read(input_alig, 'fasta')

alignment=list(align)

first_alignment=split(str(alignment[0].seq))
second_alignment=split(str(alignment[1].seq))

# Convert sequences to lists
data={'Ch_v2_0_2':first_alignment,'EuSprat':second_alignment}

# Convert data to a data.frame
df=pd.DataFrame(data)

# Remove all positions with gaps in the Atlantic herring
#https://stackoverflow.com/questions/28679930/how-to-drop-rows-from-pandas-data-frame-that-contains-a-particular-string-in-a-p
df2=df[~df.Ch_v2_0_2.str.contains("-")]

df2 = df2.copy()
# Convert gaps to missing data
df2.loc[:, ('EuSprat')] = df2.loc[:, ('EuSprat')].str.replace('-','N')

# Convert any letters to upper case
df2.loc[:, ('EuSprat')] = df2.loc[:, ('EuSprat')].str.upper()
df2.loc[:, ('Ch_v2_0_2')] = df2.loc[:, ('Ch_v2_0_2')].str.upper()

# I think I am going to create a geno file for now:

lines=input_info.readlines()
line=lines[0].split()

chromosome=[line[1]]*(df2.shape[0])
df2.loc[:, ('#CHROM')]=chromosome

positions=list(range(int(line[2]),int(line[3])+1))
df2.loc[:, ('POS')]=positions

# Change order of columns in dataframe

df2 = df2[['#CHROM','POS','Ch_v2_0_2','EuSprat']]

df2.to_csv(output,sep='\t',index=False)

input_alig.close()
input_info.close()
output.close()
