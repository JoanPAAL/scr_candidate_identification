#!/usr/bin/env python 

"""Extends the coordinates on a gff3 file to the next STOP codon downstream.
   If there is no STOP codon it returns the full sequence (default_full = True)     
"""

# This script should be placed in drac: 
# $ scp RIBOSEQ/PYTHON_SCRIPTS/3_utr_generator.py jpallares@161.116.70.109:SCRIPTS/ 
# Usage: 3_utr_generator.py file1.gff3 file2.fa
# Example: time SCRIPTS/3_utr_generator.py DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/3_nt_extended_unpacked_subset_culex_tarsalis.gff3 DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/Culex-tarsalis_knwr_CONTIGS_CtarK1.fa


# Imports #
import pyranges as pr
from pyfaidx import Fasta
import extend_orfs as eo
import pandas as pd
import sys

# Execution #

## Load gff3 and Fasta files
arg = sys.argv
g = arg[1]
fasta_path = arg[2]

#fs = Fasta(fasta_path)
p_df = pr.read_gff3(g, as_df=True) #gff3 to Dataframe 
#print('This is initial p_df:\n',p_df)


## Generate an Output in gff3 format
l = arg[1].split("/")
path = ''

# Obtain the path without file name:
for element in map(lambda x: x+"/", l[:-1]): # Maybe this is 
    path += element                          # unnecessarily complex

ig = l[-1] # Input Gff3
og = path + '3_utr_' + ig # Output gff3


## Extend up to the next ORF or until the end of the sequence:
# No need to provide cds_id; 
# after 3_nt_extend.py there is already a single instance per ID
eo.extend_orfs(p=p_df, as_df=True, fasta_path=fasta_path, default_full=True,  
               direction='down', o=og)
