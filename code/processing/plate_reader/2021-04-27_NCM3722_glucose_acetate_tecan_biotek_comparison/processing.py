#%%
import numpy as np 
import pandas as pd 

DATA_PATH = '../../../../data/plate_reader/2021-04-27_NCM3722_glucose_acetate_tecan_biotek_comparison/'

# Define the wells
glucose_wells = [f'{ROW}{COL}' for COL in np.arange(3, 11) for ROW in ['C', 'E']]
acetate_wells = [f'{ROW}{COL}' for COL in np.arange(3, 11) for ROW in ['D', 'F']]

# Load the tecan data 
tecan = pd.read_csv(f'{DATA_PATH}/raw/2021-04-27_NCM3722_glucose_acetate_tecan.csv')

# Load and annotate the bmg data
biotek = pd.read_csv(f'{DATA_PATH}/raw/2021-04-27_NCM3722_glucose_acetate_biotek.csv')

# Format the tecan
melted = tecan.melt('Cycle Nr.')

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
tecan_tidy = pd.concat(dfs, sort=False)
tecan_tidy.dropna(inplace=True)
tecan_tidy['instrument'] = 'Tecan Spark'

#%%
# Format the biotek
melted = biotek.melt(['Time', 'T'])
melted.dropna(inplace=True)

# Rename
melted.rename(columns={'Time':'time_s', 'T':'temp_C', 'variable':'well', 'value':'od_600nm'},
            inplace=True)
# Convert time to datetime
melted['time_s'] = pd.to_timedelta(melted['time_s'])
melted['time_s'] = melted['time_s'].dt.total_seconds()
melted['time_s'] -= melted['time_s'].min()
melted['instrument'] = 'Biotek Epoch2'

# Concatenate
tidy = pd.concat([tecan_tidy, melted], sort=False)


#%%
# Annotate wells
tidy['medium'] = 'water'
tidy['sample'] = 'blank'
tidy['replicate'] = 0
for i, w  in enumerate(glucose_wells):
    tidy.loc[tidy['well']==w, 'replicate'] = i + 1
    tidy.loc[tidy['well']==w, 'sample'] = 'NCM3722'
    tidy.loc[tidy['well']==w, 'medium']  = 'glucose'

for i, w  in enumerate(acetate_wells):
    tidy.loc[tidy['well']==w, 'replicate'] = i + 1
    tidy.loc[tidy['well']==w, 'sample'] = 'NCM3722'
    tidy.loc[tidy['well']==w, 'medium']  = 'acetate'

#%%
tidy.to_csv(f'{DATA_PATH}/processed/2021-04-27_NCM3722_glucose_acetate_tecan_biotek_comparison.csv', index=False)
# %%
