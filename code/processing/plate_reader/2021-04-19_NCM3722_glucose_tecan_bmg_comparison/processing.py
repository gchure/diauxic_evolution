#%%
import numpy as np 
import pandas as pd 

DATA_PATH = '../../../../data/plate_reader/2021-04-19_NCM3722_glucose_tecan_bmg_comparison/'

# Define the wells
tecan_wells = [f'{ROW}{COL}' for COL in np.arange(3, 11) for ROW in ['D', 'E']]
bmg_wells = [f'{ROW}{COL:02d}' for COL in np.arange(3, 11) for ROW in ['D', 'E']]

# Load the tecan data 
tecan = pd.read_csv(f'{DATA_PATH}/raw/2021-04-19_NCM3722_glucose_tecan_spark.csv')

# Load and annotate the bmg data
bmg = pd.read_csv(f'{DATA_PATH}/raw/2021-04-19_NCM3722_glucose_bmg_omega.csv',
                header=None)

cols = np.arange(len(bmg.values[0])) + 1
cols = [str(c) for c in cols]
bmg.columns = cols
labels = ['Time [s]', 'Temp. [°C]']
cols = []
for i, n in enumerate(bmg['1'].values):
    split = n.split(':')
    if len(split) != 1:
        labels.append(split[0])
        cols.append(split[1])
    else:
        cols.append(n)

bmg['Cycle Nr.'] = labels
bmg['1'] = cols

# Tidy the datasets
tidy_dfs = []
for name, inst, wells in zip(['Tecan Spark', 'BMG Omega'], [tecan, bmg],
                             [tecan_wells, bmg_wells]):
    melted = inst.melt('Cycle Nr.')

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
    tidy.dropna(inplace=True)
    tidy['instrument'] = name

    # Annotate wells
    tidy['sample'] = 'blank'
    tidy['replicate'] = 0
    for i, w  in enumerate(wells):
        tidy.loc[tidy['well']==w, 'replicate'] = i + 1
        tidy.loc[tidy['well']==w, 'sample'] = 'NCM3722'
    tidy_dfs.append(tidy)

merged = pd.concat(tidy_dfs, sort=False)
merged.to_csv(f'{DATA_PATH}/processed/2021-04-19_NCM3722_glucose_tecan_bmg_comparison_tidy.csv', index=False)
# %%
