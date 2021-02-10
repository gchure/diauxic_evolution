#%%
import numpy as np 
import pandas as pd
import tqdm

# Load the read data and tidy
cfu = pd.read_csv('../../../data/Favate2021/raw/TableS6_molecules_per_cfu.csv')

cfu.loc[cfu['repl']=='rep1', 'replicate'] = 1
cfu.loc[cfu['repl']=='rep2', 'replicate'] = 2
cfu['replicate'] = cfu['replicate'].astype(int)

# Look only at the transcriptional stuff
cfu = cfu[cfu['seqtype']=='rna']

# Duplicate the column for override of accession numbers
cfu['gene_name'] = cfu['target_id']

# Load the ecoli gene list to map accession number to gene name
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_index.csv')

gene_names = []
cog_class = []
cog_letter = []
cog_desc = []
for g in tqdm.tqdm(cfu['gene_name'].values):
    if 'ECB_' in g:
        sel = colicogs[colicogs['accession_number'].str.lower() ==  g.lower()]
    else:
        sel = colicogs[colicogs['name'].str.lower() == g.lower()]   
    if len(sel) > 0:        
        gene_names.append(sel['common_name'].values[0])
        cog_class.append(sel['cog_class'].values[0])
        cog_letter.append(sel['cog_letter'].values[0])
        cog_desc.append(sel['cog_desc'].values[0])
    else: 
        gene_names.append(g)
        cog_class.append('Not Assigned')
        cog_letter.append('Not Assigned')
        cog_desc.append('Not Assigned')   

cfu['gene_name'] = gene_names
cfu['cog_class'] = cog_class
cfu['cog_letter'] = cog_letter
cfu['cog_desc'] = cog_desc


#%%
# Restrict to only useful things
cfu.to_csv('../../../data/Favate2021/processed/Favate2021_reads_tidy.csv', 
             index=False)


# %%
