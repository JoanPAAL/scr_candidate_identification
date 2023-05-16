#!/usr/bin/env python 

"""  
"""

# This script should be placed in drac: 
# $ scp RIBOSEQ/PYTHON_SCRIPTS/candidate_filter.py jpallares@161.116.70.109:SCRIPTS/ 
# Usage: candidate_filter.py feature_counts.out fpkm.tsv fpkm_drop_off.tsv read_distribution.tsv second_orfs.gff3
# Example: time SCRIPTS/candidate_filter.py DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/featurecounts_merged.out DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/fpkm.tsv DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/fpkm_drop_off.tsv DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/read_distribution_featurecounts_merged.tsv DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/only_second_orf.gff3
# Another example: time SCRIPTS/candidate_filter.py DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/subset_featurecounts_merged.out  DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/fpkm.tsv DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/fpkm_drop_off.tsv DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/read_distribution_featurecounts_merged.tsv DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/only_second_orf.gff3

# Imports #
import pandas as pd
import sys
import pyranges as pr
import extend_orfs as eo

# Definition # 

def table_reader(table):
    """Processes the featureCounts output 
       into a convenient Pandas DF.
    """

    s = ''
    c = 0
    with open(table, 'r') as f:
        for line in f:
            line = line.split("\t")

            if line[0] == 'Geneid':
                line = [i.split("/")[0] for i in line]
                df = pd.DataFrame(columns=line)

            elif line[0][0] != '#':
                line = [i.strip() for i in line]
                line[1] = line[1].split(";")[0] # To format Chr
                line[4] = line[4].split(";")[0] # To format Strand
                df.loc[c] = line
                c += 1

    df.drop(columns=['Start','End'],inplace=True) # I think they are not needed
    df["Length"] = pd.to_numeric(df["Length"])

    return df


def transcript_filter(df, fpkm, fpkm_drop, rd):
    """Filters the transcripts in each DataFrame and returns 
       SCR candidates fulfilling all criteria for all of them:
           1/ Length
           2/ Expression in at least one sample
           3/ Main ORF - Second ORF Ratio
           4/ Drop Off ratio 
           5/ Read Length Distribution

    """
#################################################################################################### NADYA SCRIPT IN R
# summarise(selection_1  = any(RT_rate > 0.001 & RT_rate < 0.5 & drop_off < 0.25 & expression > 0.2 & length > 15 & first_quartile_of_the_length >= index_of_first_nonzero & index_of_last_nonzero >= last_quartile_of_the_length &  percentage_of_non_zero_bp >0.1)) 
 
    t_df = df[df['Geneid'].str.contains("_second_orf")]# Temporary DataFrame
    l = list(df.columns)

    t_df = t_df.groupby("Geneid").apply(lambda x: samples_to_rows(x))
    t_df.reset_index(inplace=True,drop=True)
    t_df["SCR_Candidate"] = 'NO'


    # 1/ Length   
    t_df.loc[t_df['Length'] >= 15, 'SCR_Candidate'] = 'YES' 

    # 2/ Expression in at least one sample
    t_df['FPKM_Second_ORF'] = t_df['Sample'].apply(lambda x: fpkm[fpkm['Geneid']
                                               .str.contains("_second_orf")][x])

    # Those samples without expression 0 will not be considered:
  #  t_df.loc[(t_df['SCR_Candidate'] == 'YES') & (t_df['FPKM_Second_ORF'] == 0), ################################ IT SHOULD BE ABOVE 0.2
    t_df.loc[(t_df['SCR_Candidate'] == 'YES') & (t_df['FPKM_Second_ORF'] < 0.2),
             'SCR_Candidate'] = 'NO'


    # 3/ Second ORF - Main ORF Rate 
    t_df['FPKM_Main_ORF'] = t_df['Sample'].apply(lambda x: fpkm[fpkm['Geneid']
                                               .str.contains("_main_orf")][x])

    # The obscure part that follows only performs this computation: ################################ IT SHOULD BE ABOVE 0.2
    #                     FPKM in Second ORF / FPKM in Main ORF
    # It also ensures FPKM in Main ORF is not 0
    t_df['Second_Main_Rate'] = t_df['Sample'].apply(lambda x: 
                   float(t_df[t_df['Sample'] == x]['FPKM_Second_ORF'].values) / 
                   float(t_df[t_df['Sample'] == x]['FPKM_Main_ORF'].values)
                   if float(t_df[t_df['Sample'] == x]['FPKM_Main_ORF'].values) 
                   > 0 else "NA")

    t_df = t_df.groupby("Sample").apply(lambda x: second_main_rate_selector(x))


    # 4/ Drop Off Rate
    t_df = t_df.groupby("Sample").apply(lambda x: drop_off_rate_selector(x, 
                                                                     fpkm_drop))


    # 5/ Read Length Distribution 
    t_df = t_df.groupby("Sample").apply(lambda x: rld_selector(x, rd))


    # Format Output
    t_df = t_df[['ID','Sample','Chr','Strand','Length','FPKM_Main_ORF', 
                 'FPKM_Second_ORF','FPKM_Drop_Off_Region','Second_Main_Rate',
                 'Drop_Off_Rate','Percentage_of_nt_Covered','First_Quartile',
                 'First_non_zero_nt','Last_Quartile','Last_non_zero_nt',
                 'Read_Distribution','SCR_Candidate']] 

    return t_df 


def samples_to_rows(df):
    """
    """

    l = df.columns

    # List of Columns (Excluding Samples)
    l_c = [i for i in l if i == 'Geneid' or i == 'Chr' or 
           i == 'Strand' or i == 'Length' or i == 'ID'] 

    # List of Samples
    l_s = [i for i in l if i != 'Geneid' and i != 'Chr' and 
         i != 'Strand' and i != 'Length' and i != 'ID'] 

    out_df = pd.DataFrame(columns=l_c)

    for sample in l_s:
        t_df = df[l_c]
        t_df['Sample'] = sample
        out_df = pd.concat([out_df, t_df])

    return out_df 
    

def second_main_rate_selector(df):
    """
    """

    candidate_status = set(df['SCR_Candidate'].values).pop()
    smr = set(df['Second_Main_Rate'].values).pop() # Second / Main Rate

    if candidate_status == 'NO':
        return df

    elif smr == 'NA':
        df['SCR_Candidate'] = 'NO'

    elif smr < 0.001 or smr > 0.5:
        df['SCR_Candidate'] = 'NO'        

    return df


def drop_off_rate_selector(df, fpkm_drop):     
    """
    """

    sample = set(df['Sample'].values).pop()
    candidate_status = set(df['SCR_Candidate'].values).pop()

    if len(fpkm_drop) == 0: # In case there is no Drop Off region available
        df['FPKM_Drop_Off_Region'] = 'NA'
        dor = 'NA' # Drop Off Rate

    else:
        value_drop = float(fpkm_drop[sample])
        value_second = float(df['FPKM_Second_ORF'])
        df['FPKM_Drop_Off_Region'] = value_drop

        if value_second > 0:
            dor = value_drop / value_second 
        else: 
            dor = 'NA'
    
    df['Drop_Off_Rate'] = dor

    if dor == 'NA' or dor > 0.25:
        df['SCR_Candidate'] = 'NO'

    return df

def rld_selector(df, rd):
    """Read Length Distribution Selector
    """

    sample = set(df['Sample'].values).pop()
    candidate_status = set(df['SCR_Candidate'].values).pop()
    rd = rd[rd['Geneid'].isin(df['Geneid'].values)]


    percentage = sample + "_Percentage_of_covered_nt"
    percentage = float(rd[percentage].values)

    first_quartile = sample + "_First_Quartile"
    first_quartile = float(rd[first_quartile].values)

    last_quartile = sample + "_Last_Quartile"
    last_quartile = float(rd[last_quartile].values)

    first_non_zero = sample + "_Index_First_Non-Zero"
    first_non_zero = float(rd[first_non_zero].values)

    last_non_zero = sample + "_Index_Last_Non-Zero"
    last_non_zero = float(rd[last_non_zero].values)
   
    df['Percentage_of_nt_Covered'] = percentage
    df['First_Quartile'] = first_quartile
    df['First_non_zero_nt'] = first_non_zero
    df['Last_Quartile'] = last_quartile
    df['Last_non_zero_nt'] = last_non_zero

    if percentage == 0:
        df['Read_Distribution'] = 'FAIL'
    elif (percentage > 10 and first_non_zero < first_quartile ############################################################ I changed this, it was set too low before
          and last_non_zero > last_quartile):
        df['Read_Distribution'] = 'OK'
    else:
        df['Read_Distribution'] = 'FAIL'

    if (candidate_status == 'YES' and 
    set(df['Read_Distribution'].values).pop() == 'FAIL'):
        df['SCR_Candidate'] = 'NO'

    return df

def summarize(df):
    """
    """

    # Select Columns
    df = df[['ID','Sample','Chr','Strand','Length','FPKM_Main_ORF',
             'FPKM_Second_ORF','FPKM_Drop_Off_Region', 
             'Percentage_of_nt_Covered','SCR_Candidate']] 

    # Process DF
    if len(df[df['SCR_Candidate'] == 'YES']) > 0:
        df = df[df['SCR_Candidate'] == 'YES'] # Just keep candidate samples
    
    df['Candidate_Samples'] = len(df)
    df = df.max().to_frame().T # Just keep the maximum values of all Samples
    df.drop(columns='Sample',inplace=True) 

    return df


# Execution #

## Load the featureCounts, FPKM, Read Distribution and gff3 tables
arg = sys.argv

df = table_reader(arg[1]) # featureCounts
df['ID'] = df['Geneid'].str.split(pat="_", expand=True)[0]

fpkm = pd.read_csv(arg[2], sep='\t')
l = list(fpkm.columns)
d = {i:(i.split(".")[0] if i != 'Geneid' else 'Geneid') for i in l}
fpkm.rename(columns=d, inplace=True) # FPKM with column names modified
fpkm['ID'] = fpkm['Geneid'].str.split(pat="_", expand=True)[0]

fpkm_drop = pd.read_csv(arg[3], sep='\t')
l = list(fpkm_drop.columns)
d = {i:(i.split(".")[0] if i != 'Geneid' else 'Geneid') 
     for i in l}
fpkm_drop.rename(columns=d, inplace=True) # FPKM Drop with column names modified
fpkm_drop['ID'] = fpkm_drop['Geneid'].str.split(pat="_", expand=True)[0]

rd = pd.read_csv(arg[4], sep='\t') # Read Distributions
rd.rename(columns={'ID':'Geneid'}, inplace=True)
rd['ID'] = rd['Geneid'].str.split(pat="_", expand=True)[0]

gff3 = pr.read_gff3(arg[5], as_df=True) # Gff3 of the Extensions


## Apply Filters
out_df = (df.groupby("ID")
          .apply(lambda x: 
                 transcript_filter(df[df['ID'].isin(x['ID'])], 
                                   fpkm[fpkm['ID'].isin(x['ID'])], 
                                   fpkm_drop[fpkm_drop['ID'].isin(x['ID'])], 
                                   rd[rd['ID'].isin(x['ID'])])))
out_df.reset_index(inplace=True,drop=True)


## Generate a Summary DataFrame
s_df = out_df.groupby("ID").apply(lambda x: summarize(x)) 
s_df.reset_index(inplace=True,drop=True)

candidate_df = s_df[s_df['SCR_Candidate'] == 'YES'] # Subselect Candidates


## Generate Output 
# Return Full Output as .tsv
l = arg[1].split("/")
path = ''

# Obtain the path without file name:
for element in map(lambda x: x+"/", l[:-1]): # Maybe this is 
    path += element                          # unnecessarily complex

fc = l[-1].strip('.out') #Feature Counts file
out = path + 'full_info_orf_candidates_' + fc + '.tsv'
out_df.to_csv(out, sep='\t', index=False) # Export output file to .tsv

# Return Summary Output as .tsv
out = path + 'summary_orf_candidates_' + fc + '.tsv' 
s_df.to_csv(out, sep='\t', index=False) # Export output file to .tsv

# Return output as gff3
# (Only Candidates Selected)
gff3 = gff3[gff3['ID'].isin(candidate_df['ID'] + "_second_orf")] 
gff3.reset_index(inplace=True,drop=True)
out = path + 'orf_candidates_' + fc + '.gff3'
eo.mod_to_gff3(gff3, path=out) # Export output file to .gff3




 
