#%%
import numpy as np 
import pandas as pd 

DATA_PATH = '../../../../data/plate_reader/2021-04-16_NCM3722_glucose_titration_tecan_comparison/'

inf = pd.read_csv(f'{DATA_PATH}/raw/2021-04-16_NCM3722_glucose_titration_tecan_infinite.csv')
spark = pd.read_csv(f'{DATA_PATH}/raw/2021-04-16_NCM3722_glucose_titration_tecan_spark.csv')

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
df['gluc_conc_mM'] = 0
df['replicate'] = 0
#%%
wells = [f'{ROW}{COL}' for COL in np.arange(3, 11) for ROW in ['C', 'D', 'E', 'F']]
concs = {n:c for n, c in zip(np.arange(3, 11), [10, 8, 6, 5, 4, 3, 2, 1])}
reps = {c:n for n, c in zip(np.arange(1, 5), ['C', 'D', 'E', 'F'])}

#%%
for w in wells:
    _rep = reps[w[0]]
    _conc = concs[int(w[1:])]
    df.loc[df['well']==w, 'replicate'] = _rep 
    df.loc[df['well']==w, 'gluc_conc_mM'] = _conc
    df.loc[df['well']==w, 'strain'] = 'NCM3722'

df.to_csv(f'{DATA_PATH}/processed/2021-04-16_NCM3722_glucose_titration_tecan_comparison_tidy.csv', 
           index=False)


# %%
