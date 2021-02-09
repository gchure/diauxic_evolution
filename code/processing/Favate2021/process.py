#%%
import numpy as np 
import pandas as pd
import diaux.io 

# Load the read data and tidy
reads = pd.read_csv('../../../data/Favate2021/raw/TableS1_read_counts.csv')
reads.loc[reads['repl']=='rep1', 'replicate'] = 1
reads.loc[reads['repl']=='rep2', 'replicate'] = 2
reads['replicate'] = reads['replicate'].astype(int)
reads.rename(columns={'target_id':'gene_name'}, inplace=True)

# Look only at the transcriptional stuff
reads = reads[reads['seqtype']=='rna']

# Restrict to only useful things
reads = reads[['line', 'replicate', 'gene_name', 'tpm']]
            
names, cogs = diaux.io.standardize_genes(reads['gene_name'].values)



# %%
