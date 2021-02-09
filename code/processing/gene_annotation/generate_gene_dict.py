#%%
import pickle
import tqdm
import pandas as pd

# Load the master gene list
colicogs = pd.read_csv('../../../data/escherichia_coli_gene_index.csv')

# Load the gene classification scheme from Balakrishnan 2021
classifier = pd.read_csv('../../../data/Balakrishnan2021/raw/genes_classification_all.csv')
classifier.loc[classifier['COG_function']=='notdefined', 'COG_function'] = 'Not Assigned'

# Instantiate the dictionary to save
colidict = {}

# Populate the dictionary using synonyms
for g, d in tqdm.tqdm(colicogs.groupby(['name']), desc='Synonyms...'):
    _colidict = {'common_name': d['common_name'].values[0],
                 'cog_letter': d['cog_letter'].values[0],
                 'cog_class': d['cog_class'].values[0],
                 'sector': ''}
    colidict[g.lower()] = _colidict

# # Populate the dictionary using the b_number
# for g, d in tqdm.tqdm(colicogs.groupby(['b_number']), desc='Accession numbers...'):
#     _colidict = {'common_name': d['common_name'].values[0],
#                  'cog_letter': d['cog_letter'].values[0],
#                  'cog_class': d['cog_class'].values[0]}
#     colidict[g.lower()] = _colidict
    
# Populate the dictionary using the common name in case it is not included as a synonym
for g, d in tqdm.tqdm(colicogs.groupby(['common_name']), desc='Common names...'):
    _colidict = {'common_name': d['common_name'].values[0],
                 'cog_letter': d['cog_letter'].values[0],
                 'cog_class': d['cog_class'].values[0],
                 'sector': ''}
    colidict[g.lower()] = _colidict


# # Given this dictionary, assign the sectors as defined in the Balakrishnan daata
# classifier['gene'] = [''.join(g.split('-')) for g in classifier['gene'].values]
# for [gene, sector], d in tqdm.tqdm(classifier.groupby(['gene', 'sector']), 
#                                    desc='Assigning sectors'):
#     # Deal with edge cases
#     if gene.lower() == 'arpB1':
#         gene = 'arpb'
#     elif gene.lower() =='insn1':
#         gene = 'insn'
#     elif gene.lower() == 'inso1':
#         continue
#     elif gene.lower() == 'rdoa':
#         continue
#     else:
#         gene = gene.split('_')[0].lower()
#     colidict[gene]['sector'] =  sector.split('-')[0].upper()

# %%
# Save the dictionary to disk 
with open('../../../diaux/package_data/coli_gene_dict.pkl', 'wb') as file:
    pickle.dump(colidict, file)

# %%
