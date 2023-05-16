#!/usr/bin/env python 

"""  
"""

# This script should be placed in drac: 
# $ scp RIBOSEQ/PYTHON_SCRIPTS/transcript_selector.py jpallares@161.116.70.109:SCRIPTS/ 
# Usage: transcript_selector.py file1.gff3 
# Example: time SCRIPTS/transcript_selector.py DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/main_orfs_found_graph_genes_regular_cluster_unpacked_subset_culex_tarsalis.gff3 DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/sequence_information_main_orfs_found_graph_genes_regular_cluster_unpacked_subset_culex_tarsalis.tsv

# Also: time SCRIPTS/transcript_selector.py DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/just_a_test.gff3 DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/sequence_information_main_orfs_found_graph_genes_regular_cluster_unpacked_subset_culex_tarsalis.tsv

# Imports #
import pyranges as pr
#from pyfaidx import Fasta
import extend_orfs as eo
import pandas as pd
import sys

# Definitions #

def transcript_counter(df):
    """
    """
    nt = len(set(df['ID'].values)) # Number of Transcripts
    df['Transcripts_per_Gene'] = nt
    return df[['Gene','Transcripts_per_Gene']].iloc[0].to_frame().T


def transcript_selector(df):
    """
    """
    if len(set(df['ID'].values)) > 1:
        df = three_prime_selector(df)
        if len(set(df['ID'].values)) > 1:
            df = length_selector(df)
            if len(set(df['ID'].values)) > 1:
                df = final_selector(df)

    return df


def three_prime_selector(df):
    """
    """
    if set(df['Strand'].values).pop() == '+': # if '+' in set()
        three_p = df['End'].max() # Most Three Prime position on + Strand
        t_df = df[df['End'] == three_p] # Temporary DataFrame 
    else: 
        three_p = df['Start'].min() # Most Three Prime position on - Strand
        t_df = df[df['Start'] == three_p] # Temporary DataFrame
    
    ids = t_df['ID']
    return df[df['ID'].isin(ids)]


def length_selector(df):
    """
    """
    t_df = (df.groupby("ID").apply(lambda x: sum(x["End"] - x["Start"])) 
                           .to_frame(name='Length'))
    t_df.reset_index(inplace=True)
    t_df = t_df[t_df['Length'] == t_df['Length'].max()]
    ids = t_df['ID']

    return df[df['ID'].isin(ids)]
 

def final_selector(df):
    """
    For those Genes that still have multiple Transcripts.
    """
    
    t_df = (df.groupby("ID").apply(lambda x: '_'.join(list(x['Start']
                                                .astype(str) + '_' + 
                                                x['End'].astype(str))))
                                   .to_frame(name='Aggregate_Positions'))
    counter = t_df['Aggregate_Positions'].value_counts()
    if len(counter) != 1: # Check if all entries are repeated
        print("We are working on gene:",set(df['Gene'].values).pop())
        print('There are multiple transcripts meeting the criteria. The first'
              ' one will be kept.') # This Warning is only shown if different
                                    # transcripts meet all criteria

    identity = df.iloc[0]['ID'] # Be it that there are multiple transcripts with 
                                # exactly the same exons or not, choose the 
                                # first one as all fulfill the criteria
    return df[df['ID'] == identity] 
                      
# Execution #

## Load gff3, fasta and tsv files
arg = sys.argv
p_df = pr.read_gff3(arg[1], as_df=True) # gff3
ti = pd.read_csv(arg[2], sep="\t") # Transcript Information


## Subset for those that have a bona fide Stop
stops = ti[ti["Stop_Codon"] == 1]
stops = stops["ID"]
p_df = p_df[p_df["ID"].isin(stops)]


## How many genes have more than one transcript?
# Generate a DataFrame with the Number of Transcripts per gene
nt_df = p_df.groupby("Gene").apply(lambda x: transcript_counter(x)) 
nt_df.reset_index(inplace=True,drop=True)
mt = len(nt_df[nt_df['Transcripts_per_Gene'] != 1]) # Multiple Transcripts
print('\nThere are ',mt,' genes with multiple transcripts.\n')


## Select a Transcript per Gene
p_df = p_df.groupby("Gene").apply(lambda x: transcript_selector(x))
p_df.reset_index(inplace=True,drop=True)


## How many genes have more than one transcript?
# Generate a DataFrame with the Number of Transcripts per gene
nt_df = p_df.groupby("Gene").apply(lambda x: transcript_counter(x)) 
nt_df.reset_index(inplace=True,drop=True)
mt = len(nt_df[nt_df['Transcripts_per_Gene'] != 1]) # Multiple Transcripts
print('\nThere are ',mt,' genes with multiple transcripts.')


## Generate an Output in gff3 format 
l = arg[1].split("/")
path = ''

# Obtain the path without file name:
for element in map(lambda x: x+"/", l[:-1]): # Maybe this is 
    path += element                          # unnecessarily complex

ig = l[-1] # Input Gff3
og = path + 'one_transcript_' + ig # Output gff3
eo.mod_to_gff3(p_df, path=og) 
