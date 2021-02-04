#%%
import pandas as pd

# Load the raw 2020 json file.
data = pd.read_json('./raw/ecoli_coglist_2020.json')

# Convert each 'results' seciton to a dataframe. 
cog_df = pd.DataFrame([])
for i in range(len(data)):
    results = data.iloc[i]['results']
    if len(results['cog']['funcats']) > 1:
        print(f"warning, {results['gene_tag']} has more than one COG")
    cog_df = cog_df.append({
        'b_number': results['gene_tag'],
        'gene_name': results['cog']['genes'][0],
        'cog_desc': results['cog']['name'],
        'cog_letter': results['cog']['funcats'][0]['name'],
        'cog_class': results['cog']['funcats'][0]['description']
    }, ignore_index=True)

cog_df

# %%
