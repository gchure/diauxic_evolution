#%%
import numpy as np
import pandas as pd

# Load the data
data = pd.read_csv('./input/2021-02-11_REL606_NCM3722_diauxie.csv')

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
for n in range(4, 10):
    if n <= 6:
        medium = '10 mM glucose + 30 mM acetate'
    else:
        medium = '0.61 mM glucose + 30 mM acetate'
    for letter, strain in zip(['D', 'E'], ['NCM3722', 'REL606']):
        tidy.loc[tidy['well'] == f'{letter}{n}', 'strain'] = strain
        tidy.loc[tidy['well'] == f'{letter}{n}', 'medium'] =  medium

# Add replicate information.
for g, d in tidy.groupby(['strain', 'medium']):
    mapper = {w:r + 1 for r, w in enumerate(d['well'].unique())}
    for k, v in mapper.items():
        tidy.loc[tidy['well']==k, 'replicate'] = v
tidy['replicate'] = tidy['replicate'].astype(int)

# Save the tidy dataframe to disk for further processing
tidy.to_csv('./output/2021-02-11_NCM_REL_diauxie_tidy.csv', index=False)


# %%

# %%
