#%%
import pandas as pd
import numpy as np
import tqdm

# Load the ecocyc gene list
ecocyc = pd.read_csv('../../../data/ecocyc/escherichia_coli_mg1655_ecocyc_annotations.txt',
                     delimiter='\t')
# Load the COGs
cog = pd.read_csv('../../../data/COG_annotations/escherichia_coli_cogs_2020.csv')


#%%
# make a look-up table of gene name and b number to common name. 
coli_gene_lut = pd.DataFrame([])
for g, d in tqdm.tqdm(ecocyc.groupby(['gene_name']), desc='processing EcoCyc gene list'):
    b_number = d['b_number'].values[0]

    # Get the ecocyc data
    common_name = d['common_name'].values[0]
    if type(common_name) == float:
        continue   

    # Fix weirdness in including italics
    common_name = common_name.replace('<i>', '')
    common_name = common_name.replace('</i>', '')
    if d['summary'].isnull().all():
        summary = 'No EcoCyc summary available'
    else:
        summary = d['summary'].values[0]

    # Get the cog data
    cog_entry = cog[cog['b_number']==b_number]
    if len(cog_entry) == 0:
        cog_class = 'Not Assigned'
        cog_letter = "Not Assigned"
        cog_desc = 'Not Assigned'
    else:
        cog_class = cog_entry['cog_class'].values[0]
        cog_letter = cog_entry['cog_letter'].values[0]
        cog_desc = cog_entry['cog_desc'].values[0]
    # Check if there are synonyms

    if ~d['synonyms'].isnull().all():
        syn = d['synonyms'].values[0]
        if '\\' in syn:
            split = ' \\ '
        else:
            split = ' // '
        synonyms = d['synonyms'].values[0].split(split)

        # If there are synonyms, add them to the dataframe individually
        for i, s in enumerate(synonyms):
            coli_gene_lut = coli_gene_lut.append({ 
                        'name': s,
                        'common_name': common_name,
                        'cog_class': cog_class,
                        'cog_letter': cog_letter,
                        'cog_desc': cog_desc,
                        'ecocyc_summary': summary, 
                        },
                        ignore_index=True)
    else:
        # If there are no synonyms, add the gene name directly.
        coli_gene_lut = coli_gene_lut.append({
                        'name': g,
                        'common_name': common_name,
                        'cog_class': cog_class,
                        'cog_letter': cog_letter,
                        'cog_desc': cog_desc,
                        'ecocyc_summary': summary
                        }, 
                        ignore_index=True)


#%%
coli_gene_lut.to_csv('../../../data/escherichia_coli_gene_index.csv', index=False)


# %%
