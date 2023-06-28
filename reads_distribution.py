#!/usr/bin/env python 

"""  
"""

# This script should be placed in drac: 
# $ scp RIBOSEQ/PYTHON_SCRIPTS/reads_distribution.py jpallares@161.116.70.109:SCRIPTS/ 
# Usage: reads_distribution.py feature_counts.out file.gff3 path_to_folders_containing_bam_files
# Example: time SCRIPTS/reads_distribution.py DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/featurecounts_merged.out DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/no_parent_data.gff3 DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/STAR_ALIGNED_OUTPUT/
# Another example: time SCRIPTS/reads_distribution.py DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/subset_featurecounts_merged.out DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/NO_PARENT_DATA_small_sample_with_main_and_second.gff3 DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/STAR_ALIGNED_OUTPUT/


# Imports #
import pandas as pd
import sys
import plastid as plid


# Definition # 

def sample_iterator(table, gff3, path_to_bam):

    """Allows to iterate through samples from table.
       For each sample it calculates different statistics on Read Distribution
       by calling the function read_statistics()       
    """

    df = pd.DataFrame() # Generate the Output DataFrame (Empty at this point)

    with open(table, 'r') as f: 
        for line in f:
            line = line.split("\t")
           
            if line[0] == 'Geneid':
                samples = [i.strip() for i in line if i != 'Geneid' and 
                           i != 'Chr' and i != 'Start' and i != 'End' and 
                           i != 'Strand' and i != 'Length']

                for sample in samples:
                    path_to_bam_sample = path_to_bam + sample
                    alignment = plid.BAMGenomeArray(path_to_bam_sample)
                    
                    # To trim n = nibble nts from each side
                    # Reads for which this cannot be applied will be ignored
                    alignment.set_mapping(plid.CenterMapFactory(nibble=12))
                    sample = sample.split("/")[0]
                    transcripts = plid.GFF3_TranscriptAssembler(gff3) 

                    # Generate a Temporary DataFrame with the read statistics
                    df_t = read_statistics(transcripts, sample, alignment) 

                    # Update df with df_t
                    if len(df) == 0:
                        df = df_t
                    else:
                        df = pd.merge(df, df_t, on='ID')       
                           
                break # The rest of the file is not relevant here
    return df


def read_statistics(transcripts, sample, alignment):

    """Computes Read Distribution Statistics for a specific sample
       For each sample the number of reads per nt is obtained.
       This allows to compute the % of nt covered and the index of the first and
       last nt with reads. 
       These indexes and the quartiles of the sequence length provide with an
       additional estimate of read distribution.
    """

    df = pd.DataFrame(columns=["ID", sample + "_Percentage_of_covered_nt", 
                               sample + "_First_Quartile",
                               sample + "_Last_Quartile", 
                               sample + "_Index_First_Non-Zero", 
                               sample + "_Index_Last_Non-Zero"]) 

    c = 0
    for transcript in transcripts: 
        name = transcript.get_name()
        lt = [name] # Temporary List

        counts = transcript.get_counts(alignment) # Get reads per nt

        non_zero_nt = [i for i in counts if i !=0]
        percentage = round(len(non_zero_nt) / len(counts) * 100, 2)

        if percentage == 0:
            first_quartile = 'NA'
            last_quartile = 'NA'
            index_of_first_nonzero = 'NA'
            index_of_last_nonzero = 'NA'

        else:
            first_quartile = (int(len(counts)) / 4) - 1 
            last_quartile = (int(len(counts)) / 4 * 3) - 1
            index_of_first_nonzero = [i for i in range(len(counts)) 
                                      if counts[i] != 0 ][0]
            index_of_last_nonzero = [i for i in range(len(counts)) 
                                     if counts[i] != 0 ][-1]

        lt.append(percentage)
        lt.append(first_quartile)   
        lt.append(last_quartile)
        lt.append(index_of_first_nonzero)
        lt.append(index_of_last_nonzero)

        df.loc[c] = lt

        c += 1 

    return df

# Execution #

## Load Files and Calculate Read Distribution Statistics
arg = sys.argv
df = sample_iterator(arg[1], arg[2], arg[3])


## Generate Output
l = arg[1].split("/")
path = ''

# Obtain the path without file name:
for element in map(lambda x: x+"/", l[:-1]): # Maybe this is 
    path += element                          # unnecessarily complex

fc = l[-1].strip('.out') #Feature Counts file
out = path + 'read_distribution_' + fc + '.tsv'

# Export output file to .tsv
df.to_csv(out, sep='\t', index=False)














