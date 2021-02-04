#%%
import requests
import pandas as pd
import json

# Define the NCBI COG url
URL = "https://www.ncbi.nlm.nih.gov/research/cog/api/cog/?organism=Escherichia_coli_K-12_sub_MG1655&format=json" 
content = requests.get(URL)
cog_json = json.loads(content.content)

cog_df = pd.DataFrame([])
# Set up the dataframe simplifying the COG information to only what we need
NEXT_PAGE = True
iter = 1
while NEXT_PAGE:
    print(f'Processing page {iter}...')
    for result in cog_json['results']:
        for funcat in result['cog']['funcats']: 
            if type(result['cog']['genes']) is not list:
                gene = result['cog']['genes']
            else:
                if len(result['cog']['genes']) > 1:
                    raise ValueError(f"more than one gene! {result['cog']['genes']}")
                elif len(result['cog']['genes']) == 0:
                    continue 
                else: 
                    gene = result['cog']['genes'][0]
            cog_dict = {'b_number': result['gene_tag'],
                    'gene': gene,
                    'cog_desc':result['cog']['name'],
                    'cog_class': funcat['description'],
                    'cog_letter': funcat['name']}
            # Update the dataframe
            cog_df = cog_df.append(cog_dict, ignore_index=True)

    # Load the next page
    if cog_json['next'] != None:
        iter += 1
        NEXT_PAGE = True
        URL = cog_json['next']
        content = requests.get(URL)
        cog_json = json.loads(content.content)
    else:
        NEXT_PAGE = False
print('finished!')
# %%
cog_df.to_csv('../../../data/COG_annotations/escherichia_coli_cogs_2020.csv', 
              index=False)

# %%
