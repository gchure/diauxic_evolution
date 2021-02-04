"""
Note that this script can only be run if BioCyc's pathway tools is installed 
locally and is being served via the following command:

./pathway-tools -lisp -python
"""
#%%
import pandas as pd
import tqdm

# Load the ecocyc gene list
ecocyc = pd.read_csv('../../../data/ecocyc/raw/ecoli_mg1655_gene_list.txt',
                     delimiter='\t')

# Load the COGs
cog = pd.read_json('../../../data/COG_annotations/raw/ecoli_coglist_2020.json')
cog.head()
#%%
# make a look-up table of gene name and b number to common name. 
coli_gene_lut = pd.DataFrame([])
for g, d in tqdm.tqdm(ecocyc.groupby(['b_number']), desc='processing EcoCyc gene list'):
    common_name = d['common_name'].values[0]
    # Check if there are synonyms
    if ~d['synonyms'].isnull().all():
        syn = d['synonyms'].values[0]
        if '\\' in syn:
            split = ' \\ '
        else:
            split = ' // '
        synonyms = d['synonyms'].values[0].split(split)
        for i, s in enumerate(synonyms):
            coli_gene_lut = coli_gene_lut.append({
                        'b_number': g,
                        'synonym': s,
                        'common_name': common_name},
                        ignore_index=True)
    else:

        coli_gene_lut = coli_gene_lut.append({
                        'b_number': g,
                        'synonym': d['gene_name'].values[0],
                        'common_name': common_name}, 
                        ignore_index=True)


coli_gene_lut.head()

# %%
