#%%
import numpy as np 
import pandas as pd
import tqdm
import diaux.io 
#%%
# Load the read data and tidy
reads = pd.read_csv('../../../data/Favate2021/raw/TableS1_read_counts.csv')
reads.loc[reads['repl']=='rep1', 'replicate'] = 1
reads.loc[reads['repl']=='rep2', 'replicate'] = 2
reads['replicate'] = reads['replicate'].astype(int)

# Look only at the transcriptional stuff
reads = reads[reads['seqtype']=='rna']

# Duplicate the column for override of accession numbers
reads['gene_name'] = reads['target_id']

# Load the ecoli gene list to map accession number to gene name
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_index.csv')
colicogs = colicogs[colicogs['strain']=='REL606']

for g in tqdm.tqdm(reads['gene_name'].unique()):
    if 'ECB_' in g:
        common_name = colicogs[colicogs['accession_number']==g]['common_name']
        if len(common_name) > 0:
            reads.loc[reads['gene_name']==g, 'gene_name'] = common_name.values[0]

# Restrict to only useful things
reads = reads[['line', 'replicate', 'gene_name', 'tpm']]

# Convert all names to  standard name and get COG/sector info
names = [''.join(g.split('_')) for g in reads['gene_name'].values] 
std_dict = diaux.io.standardize_genes(names)



# %%
