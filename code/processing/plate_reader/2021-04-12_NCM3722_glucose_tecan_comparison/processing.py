#%%
import numpy as np 
import pandas as pd 

DATA_PATH = '../../../../data/plate_reader/2021-04-12_NCM3722_glucose_tecan_comparison/'

inf = pd.read_csv(f'{DATA_PATH}/raw/2021-04-12_NCM3722_glucose_tecan_infinite.csv')
spark = pd.read_csv(f'{DATA_PATH}/raw/2021-04-12_NCM3722_glucose_tecan_spark.csv')

inst_dfs = [] 
for data, name in zip([inf, spark], ['infinite', 'spark']):
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
    tidy.dropna(inplace=True)
    tidy['instrument'] = name
    inst_dfs.append(tidy)

# Combine into a single huge dataframe
df = pd.concat(inst_dfs, sort=False)

# Assign the wells
df['strain'] = 'blank'
df['replicate'] = 0

reps = [f'{ROW}{COL}' for COL in np.arange(3, 9) for ROW in ['C', 'D', 'E', 'F']]
reps.append('F9')
for i, r in enumerate(reps):
    i += 1
    df.loc[df['well']==r, 'replicate'] = i
    df.loc[df['well']==r, 'strain'] = 'NCM3722'

df.to_csv(f'{DATA_PATH}/processed/2021-04-12_NCM3722_tecan_comparison_tidy.csv', 
           index=False)
# %%
