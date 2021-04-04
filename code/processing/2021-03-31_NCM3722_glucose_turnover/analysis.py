#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import scipy.stats
import altair as alt
import cmdstanpy
import arviz as az
import seaborn as sns 
sns.set()

# Load the various datasets 
DATA_PATH = '../../../data/2021-03-31_NCM3722_glucose_turnover/processed'
growth_data = pd.read_csv(f'{DATA_PATH}/2021-03-31_NCM3722_glucose_growth_data.csv')
cal_data = pd.read_csv(f'{DATA_PATH}/2021-03-31_glucose_acetate_calibration.csv')
samp_data = pd.read_csv(f'{DATA_PATH}/2021-03-31_NCM3722_glucose_turnover_peaks.csv')
# %%
# Take the growth data and link the compound peak area to sample id
turnover_dfs = []
for g, d in growth_data.groupby(['replicate', 'sample_idx']):
    if g[1] == 0:
        continue
    # Get the corresponding replicate and sample ID from the turnover data
    _peak = samp_data[(samp_data['replicate']==g[0]) & 
                      (samp_data['sample_idx']==g[1])]
    for _g, _d in _peak.groupby(['compound']):
        d[f'{_g}_area'] = _d['area'].values[0]
    turnover_dfs.append(d)
turnover_df = pd.concat(turnover_dfs, sort=False)
turnover_df.fillna(0, inplace=True)
# Standardize to the average phsophate area of the "right" peaks
avg_phos = turnover_df[turnover_df['phosphate_area'] >= 6E6]['phosphate_area'].mean()
turnover_df['correction_factor'] = turnover_df['phosphate_area'].values / avg_phos
turnover_df['phosphate_area'] /= turnover_df['correction_factor'].values
turnover_df['glucose_area'] /= turnover_df['correction_factor'].values
turnover_df['acetate_area'] /= turnover_df['correction_factor'].values
#%% 
# Look at the posphate peak area 
alt.Chart(turnover_df).mark_point().encode(
        x='od_600nm:Q',
        y='phosphate_area:Q',
        color = 'replicate:N'
)
# %%
# Do a simple linear regression to get the calibration details
params = {}
for g, d in cal_data.groupby(['compound']):
    if g  != 'phosphate':
        popt = scipy.stats.linregress(d['concentration_mM'].values, d['area'].values)
        params[g] = {'slope':popt[0],
                     'intercept':popt[1]}
# %%
# Convert the integrated areas to concentration
for carb in ['glucose', 'acetate']:
    slope = params[carb]['slope']
    inter = params[carb]['intercept']
    turnover_df[f'{carb}_conc_mM'] = (turnover_df[f'{carb}_area'].values + inter) / slope
    turnover_df.loc[turnover_df[f'{carb}_conc_mM'] < 0, f'{carb}_conc_mM'] = 0
    
# %%
from altair_saver import save
chart = alt.Chart(turnover_df).mark_point(filled=True, size=80).encode(
        x=alt.X('od_600nm:Q', title='optical density [a.u.]',
                scale=alt.Scale(zero=False)),
        y=alt.Y('acetate_conc_mM:Q', title='acetate concentration [mM]'),
        color =alt.Color('replicate:N', title='biological replicate')
)
save(chart, './output/NCM3722_acetate_production.pdf')
chart
# %%
chart = alt.Chart(turnover_df).mark_point(filled=True, size=80).encode(
        x=alt.X('glucose_conc_mM:Q', title='glucose concentratoin [mM]',
                scale=alt.Scale(zero=False)),
        y=alt.Y('acetate_conc_mM:Q', title='acetate concentration [mM]',
                scale=alt.Scale(zero=False)),
        color =alt.Color('replicate:N', title='biological replicate')
)
save(chart, './output/NCM3722_glucose_acetate_turnover.pdf')
chart
# %%
chart = alt.Chart(turnover_df).mark_point(filled=True, size=80).encode(
        x=alt.X('od_600nm:Q', title='optical density [a.u.]'),
        y=alt.Y('glucose_conc_mM:Q', title='glucose concentration [mM]',
                scale=alt.Scale(domain=[6, 11])),
        color =alt.Color('replicate:N', title='biological replicate')
)
save(chart, './output/NCM3722_glucose_consumption.pdf')
chart
# %%

# %%
