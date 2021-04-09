#%% 
sys.path.insert(0, '../../../../')
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import diaux.viz
import cremerlab.growth
colors, palette = diaux.viz.matplotlib_style()

# Define the data path and load raw data
DATA_PATH = '../../../../data/metabolite_turnover/2021-04-04_REL606_glucose_turnover'
raw_data = pd.read_csv(f'{DATA_PATH}/raw/2021-04-04_REL606_glucose_growth_curve.csv')
data, params, plot = cremerlab.growth.infer_growth_rate(raw_data, 
                                                  od_bounds=[0.039, 0.41],
                                                  groupby=['replicate'],
                                                  viz=True)
fig, ax = plot
for ext in ['pdf', 'png']:
    plt.savefig(f'./output/2021-04-04_REL_growth_curves.{ext}')
#%%
# Transforme the boolean sampled column to annotated time points
dfs = []
for g, d in data.groupby(['replicate']):
        d['time_idx'] = 0
        sample_idx = []
        iter = 0
        for s in d['sampled'].values:
                if s == True:
                        iter += 1
                        sample_idx.append(iter)
                else:
                        sample_idx.append(0)
        d['time_idx'] = sample_idx
        dfs.append(d)

data = pd.concat(dfs, sort=False)
data.drop(columns=['sampled'], inplace=True)

# Save the pruned data and growth rate parameters
data.to_csv(f'{DATA_PATH}/processed/2021-04-04_REL606_glucose_growth.csv')
params.to_csv(f'{DATA_PATH}/processed/2021-04-04_REL606_glucose_growth_parameters.csv')
#%%
