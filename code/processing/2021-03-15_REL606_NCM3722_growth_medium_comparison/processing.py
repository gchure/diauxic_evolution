#%%
import numpy as np
import pandas as pd

# Load the data
data = pd.read_csv('../../../data/2021-03-15_REL606_NCM3722_growth_medium_comparison/raw/2021-03-15_REL606_NCM3722_NC_DM_comparison.csv')

#%%
# DO some serious tidying
melted = data.melt('Cycle Nr.')

# Get the time indices
time = melted[melted['Cycle Nr.']=='Time [s]']
time.sort_values(by='variable', inplace=True)
time = time['value'].values

# Get the temperature indices
temp = melted[melted['Cycle Nr.']=='Temp. [°C]']
temp.sort_values(by='variable', inplace=True)
temp = temp['value'].values

# get the well info
dfs = []
_melted = melted[(melted['Cycle Nr.'] != 'Time [s]') & 
                 (melted['Cycle Nr.'] != 'Temp. [°C]')]
for g, d in _melted.groupby(['Cycle Nr.']):
    d.sort_values(by='variable', inplace=True)
    d['time_s'] = time
    d['temp_C'] = temp
    d.rename(columns={'Cycle Nr.': 'well', 
                      'value':'od_600nm'}, inplace=True)
    d.drop(columns=['variable'], inplace=True)
    dfs.append(d)
tidy = pd.concat(dfs, sort=False)

# Add identifier for the strain
tidy['strain'] = 'blank'
tidy['medium'] = 'blank'
tidy['glucose_mM'] = 0 
for n in range(3, 11):
    for letter, strain in zip(['C', 'D', 'E', 'F'], 
                              ['NCM3722', 'NCM3722', 'REL606', 'REL606']):
        if (letter =='C') | (letter=='E'):
            glucose = 10
        else: 
            glucose = 1.39

        if n <= 6:
            medium = 'N-C- base'
        else:
            medium = 'DM base'        

        tidy.loc[tidy['well'] == f'{letter}{n}', 'strain'] = strain
        tidy.loc[tidy['well'] == f'{letter}{n}', 'medium'] =  medium
        tidy.loc[tidy['well'] == f'{letter}{n}', 'glucose_mM'] =  glucose

# Add replicate information.
for g, d in tidy.groupby(['strain', 'medium']):
    mapper = {w:r + 1 for r, w in enumerate(d['well'].unique())}
    for k, v in mapper.items():
        tidy.loc[tidy['well']==k, 'replicate'] = v
tidy['replicate'] = tidy['replicate'].astype(int)

# Save the tidy dataframe to disk for further processing
tidy.to_csv('../../../data/2021-03-15_REL606_NCM3722_growth_medium_comparison/processed/2021-03-15_REL606_NCM3722_NC_DM_comparison_tidy.csv', index=False)


# %%

# %%
