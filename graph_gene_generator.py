#!/usr/bin/env python 

"""
    This script groups Transcripts into Gene Units according to their genomic
    overlap.
    The input is a gff3 file with information on transcript overlap 
    (i.e. Cluster).
    Upon this, graphs are generated through the NetworkX package. A graph
    consists of all transcripts that cluster at some point and is considered a
    Gene Unit.
    The output is a gff3 file with the Gene Unit information stored as Gene. 

    The essential steps are:    
                             1/ Identify all the Transcript IDs present on each
                                Cluster
                             2/ Intersect each Cluster with each other to see
                                which share Transcript IDs
                             3/ Build Graphs in which nodes are clusters and 
                                edges represent the clusters (nodes) share 
                                Transcript IDs
                             4/ Give to each Graph a unique Gene Name
"""

# This script should be placed in drac: 
# $ scp RIBOSEQ/PYTHON_SCRIPTS/graph_gene_generator.py jpallares@161.116.70.109:SCRIPTS/ 
# Usage: graph_gene_generator.py file1.gff3 
# Example: time SCRIPTS/graph_gene_generator.py DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/regular_cluster_unpacked_subset_culex_tarsalis.gff3 # This takes around 9 minutes
# Some tests: time SCRIPTS/graph_gene_generator.py DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/CULEX_TARSALIS_GENOME/shuf_prova_gene_generator.gff3

# Imports #
import pyranges as pr
import extend_orfs as eo
import pandas as pd
import sys
import networkx as nx
import itertools as it
import matplotlib.pyplot as plt

# Definition #

def id_to_gene(df, d):
    """
    """
    identity = set(df['ID'].values).pop()
    df['Gene'] = d[identity]

    return df   

# Execution #

## Load gff3 and Fasta files
arg = sys.argv
g = arg[1]
p_df = pr.read_gff3(g, as_df=True) #gff3 to Dataframe 


## Clean the DataFrame
p_df = p_df[['Chromosome', 'Start', 'End', 'Strand', 'ID', 'Parent', 
       'Feature', 'Cluster']][p_df.Feature=='CDS']
# So the Pyragnes object is correctly labelled as stranded:
p_df['Strand'] = p_df['Strand'].astype("string") 
print("This is initial p_df:\n",p_df)


## Generate Cluster-ID dictionary:
dci = {} # Cluster-ID Dictionary 

# Store cluster in a list, then iterate through p_df through it
clusters = p_df['Cluster'].unique().tolist() 

for cluster in clusters:
    t_df = p_df.loc[p_df['Cluster'] == cluster] # Temporary DataFrame
    t_s = set(t_df['ID'].values) # Temporary Set

    dci[cluster] = t_s


## Working with graphs:

# Generate graph and add nodes
g = nx.Graph()
g.add_nodes_from(dci)

# if node1 and node2 intersect in their transcripts, grow an edge between them
print('\nStarting with the edges!')
for combination in it.combinations(dci,2):
    print('This is a combiniation:\n',combination)
    intersection = dci[combination[0]].intersection(dci[combination[1]])
    print('This is an intersection:\n',intersection)
    if len(intersection) > 0:
#        g.add_edge(combination[0],combination[1], length=100000) # length is for plot
        g.add_edge(combination[0],combination[1]) 
print('Edges are done!')

# Show Componenents

components = nx.connected_components(g)
dig = {} # ID-Gene Dictionary
cm = 1 # Multiple transcript Counter
cs = 1 # Single transcript Counter

print('Starting to build dig.')
for component in components:
    t_l = [dci[i] for i in component] # Temporary List
    t_s = set(it.chain.from_iterable(t_l)) # Temporary Set
    print("This is a component:\n",component)
    print("This is t_s:\n",t_s)
    if len(t_s) == 1:
        gene_name = "s_" + str(cs)
        cs += 1 # Increase Gene Number for the next one
    else:
        gene_name = "m_" + str(cm)
        cm += 1 # Increase Gene Number for the next one
    for identity in t_s: 
        dig[identity] = gene_name 

print('dig is built.')
#print("This is dig:\n",dig)

# Graph Visualization
#plt.figure(figsize=(16, 10))
#pos = nx.spring_layout(g, k=0.3, iterations=20) # k is for node distance
#nx.draw(g, pos, with_labels=True, node_color='orange', node_size=1500, width=3)
#plt.savefig("graph_test.png")


## Give each row its corresponding Gene name using the ID-Gene Dictionary
print('About to enter lambda.')
p_df['Gene'] = "0"
p_df = p_df.groupby(by='ID').apply(lambda x: id_to_gene(x, dig))
print("\nThis is p_df towards the end:\n",p_df)


## Return output as gff3
l = arg[1].split("/")
path = ''

# Obtain the path without file name:
for element in map(lambda x: x+"/", l[:-1]): # Maybe this is 
    path += element                          # unnecessarily complex

ig = l[-1] # Input gff3
og = path + 'graph_genes_' + ig # Output gff3
eo.mod_to_gff3(p_df, path=og) 
