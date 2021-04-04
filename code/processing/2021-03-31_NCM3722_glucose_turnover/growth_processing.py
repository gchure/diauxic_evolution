#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
# import diaux.viz 
import cremerlab.growth
# colors, palette = diaux.viz.altair_style()

DATA_PATH = '../../../data'

# Load the raw data 
data = pd.read_csv(f'{DATA_PATH}/2021-03-31_NCM3722_glucose_turnover/raw/2021-03-31_NCM3722_glucose_growth.csv')
data_df, param_df, plot = cremerlab.growth.infer_growth_rate(data, viz=True,
                                                             print_params=False,
                                                             groupby=['replicate']) 

fig, ax = plot
for ext in ['pdf', 'png']:
        plt.savefig(f'./output/2021-03-31_NCM3722_glucose_growth_curves.{ext}', 
                        bbox_inches='tight')

#%%
# Convert the 'sampled' boolean to time point
dfs = []
for g, d in data_df.groupby(['replicate']):
        d['sample_idx'] = 0
        sample_idx = []
        iter = 0
        for s in d['sampled'].values:
                if s == True:
                        iter += 1
                        sample_idx.append(iter)
                else:
                        sample_idx.append(0)
        d['sample_idx'] = sample_idx
        dfs.append(d)
data_df = pd.concat(dfs, sort=False)
data_df.to_csv(f'{DATA_PATH}/2021-03-31_NCM3722_glucose_turnover/processed/2021-03-31_NCM3722_glucose_growth_data.csv')
param_df.to_csv(f'{DATA_PATH}/2021-03-31_NCM3722_glucose_turnover/processed/2021-03-31_NCM3722_glucose_growth_params.csv')


# %%
