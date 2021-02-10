#%%
import pandas as pd
import numpy as np
import tqdm

# Load the ecocyc gene list
strains = ['MG1655', 'REL606']

#%%
cog = pd.read_csv('../../../data/COG_annotations/escherichia_coli_cogs_2020.csv')
cog['gene'] = cog['gene'].str.lower()

#%%
# Instantiate the storage dataframe
coli_gene_lut = pd.DataFrame([])

# Instantiate a look-up table for gene to MG1655 b number (where possible)
b_numbers_dict = {}
# ITerate through each strain 
for strain in strains:
    ecocyc = pd.read_csv(f'../../../data/ecocyc/escherichia_coli_{strain}_ecocyc_annotations.txt',
                         delimiter='\t')
    # Iterate through each gene names 
    for g, d in tqdm.tqdm(ecocyc.groupby(['gene_name']), 
        desc=f'Processing EcoCyc gene list for {strain}'):

        # Get identifying accession numbers and products
        accession = d['ecocyc_accession'].values[0]
        product = d['gene_product'].values[0]
        common_name = d['common_name'].values[0]
        if type(common_name) == float:
            continue   

        # Fix weirdness in including italics
        common_name = common_name.replace('<i>', '')
        common_name = common_name.replace('</i>', '')   

        # Store the B number for MG1655 and map where possible
        if strain == 'MG1655':
           b_number = d['b_number'].values[0]
           if len(d['b_number']) > 0:
              b_numbers_dict[g] = b_number
              b_numbers_dict[common_name] = b_number
           else:
               b_number = 'Not Assigned'
               b_numbers_dict[g] = 'Not Assigned'
               b_numbers_dict[common_name] = 'Not Assigned'
        if strain == 'REL606':
            try:
                b_number = b_numbers_dict[g]
            except KeyError:
                b_number = 'Not Assigned'

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
                            'strain': strain,
                            'b_number': b_number,
                            'name': s,
                            'accession_number': accession,
                            'common_name': common_name,
                            'product': product},
                            ignore_index=True)
        else: 
            # If there are no synonyms, add the gene name directly.
            coli_gene_lut = coli_gene_lut.append({
                            'strain': strain,
                            'name': g,
                            'accession_number': accession,
                            'common_name': common_name,
                            'product': product,
                            'b_number': b_number}, 
                            ignore_index=True)
        coli_gene_lut = coli_gene_lut.append({
                            'strain': strain,
                            'name': common_name,
                            'accession_number': accession,
                            'common_name': common_name,
                            'product': product,
                            'b_number': b_number}, 
                            ignore_index=True)

#%%
# Iterate through the genes and link.
for i, g in tqdm.tqdm(enumerate(coli_gene_lut['name'].values), 
                      desc='Linking COGs...'):
    # Try first linking with b_number
    d = coli_gene_lut[coli_gene_lut['name']==g]
    _b_number = d[d['b_number']!='Not Assigned']['b_number']
    COG = False
    if len(_b_number) > 0:
        cog_entry = cog[cog['b_number']==_b_number.values[0]]
        if len(cog_entry) > 0:
            cog_letter = cog_entry['cog_letter'].values[0] 
            cog_class = cog_entry['cog_class'].values[0] 
            cog_desc = cog_entry['cog_desc'].values[0] 
            COG = True
        else:
            cog_letter = 'Not Assigned'
            cog_class = 'Not Assigned'
            cog_desc = 'Not Assigned'
    else: 
        cog_entry = cog[cog['gene'] == g.lower()]
        _b_number = 'Not Assigned'
        if len(cog_entry) != 0:
            cog_letter = cog_entry['cog_letter'].values[0] 
            cog_class = cog_entry['cog_class'].values[0] 
            cog_desc = cog_entry['cog_desc'].values[0] 
            COG = True
    if COG == False:
        cog_letter = 'Not Assigned'
        cog_class = 'Not Assigned'
        cog_desc = 'Not Assigned'

    coli_gene_lut.loc[coli_gene_lut['name'] == g, 'cog_letter'] = cog_letter
    coli_gene_lut.loc[coli_gene_lut['name'] == g, 'cog_class'] = cog_letter
    coli_gene_lut.loc[coli_gene_lut['name'] == g, 'cog_desc'] = cog_letter

#%%
# Group by the cases where there is more than 1 cog designation for each common name
coli_gene_cogs = []
for g, d in tqdm.tqdm(coli_gene_lut.groupby(['common_name']), 
                      desc='Consolidating COG assignments...'):
    if len(d['cog_letter'].unique()) > 1:
        _assigned = d[d['cog_letter']!= 'Not Assigned']
        d['cog_letter'] = _assigned['cog_letter'].values[0]
        d['cog_class'] = _assigned['cog_class'].values[0]
        d['cog_desc'] = _assigned['cog_desc'].values[0]  
    coli_gene_cogs.append(d)
coli_gene_cogs = pd.concat(coli_gene_cogs, sort=False)



#%%
coli_gene_cogs.to_csv('../../../data/escherichia_coli_gene_index.csv', index=False)


# %%
