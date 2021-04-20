#%%
import numpy as np
import pandas as pd 
import altair as alt
import scipy.stats
from altair_saver import save
import diaux.viz
colors, palette = diaux.viz.altair_style()
alt.data_transformers.disable_max_rows()
# %%
DATA_PATH = '../../../../data/plate_reader/2021-04-16_NCM3722_glucose_titration_tecan_comparison/processed/'
data = pd.read_csv(f'{DATA_PATH}/2021-04-16_NCM3722_glucose_titration_tecan_comparison_tidy.csv')

# Convert time to hours 
data['time_hr'] = data['time_s'] / 3600
data.drop(columns='time_s', inplace=True)
# Compute the sub
dfs = []
for g, d in data.groupby(['instrument', 'time_hr']):
    avg_blank = d[d['strain']=='blank']['od_600nm'].mean()
    d['od_600nm_sub'] = d['od_600nm'].values -  avg_blank
    dfs.append(d)
subbed = pd.concat(dfs, sort=False)
subbed = subbed[subbed['strain']!='blank']
#%%
# Plot the averages
out = []
for instrument in ['infinite', 'spark']:
    _data = subbed[subbed['instrument']==instrument]
    avg = alt.Chart(_data).mark_line().encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('mean(od_600nm_sub):Q', title='optical density [a.u.]',
        scale=alt.Scale(type='log')),

        color=alt.Color('gluc_conc_mM:O', title='glucose [mM]'))

    err = alt.Chart(_data).mark_errorband(extent='stdev').encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]'),
        color=alt.Color('gluc_conc_mM:O', title='glucose [mM]'))


    layer = (avg + err).properties(title=f'Tecan {instrument}')
    out.append(layer)

result = out[0] | out[1]
result
# save(result, 'output/2021-04-16_glucose_titration_tecan_comparison_avg.pdf')

# %%
# Given the results above, only work with the tecan spark
grouped = subbed.groupby(['instrument', 'time_hr', 'gluc_conc_mM']).mean().reset_index()
grouped = grouped[grouped['instrument']=='spark']


# Find the first point where the first derivative is about zero
fit_dfs, crop_dfs= [], []
yield_df = pd.DataFrame([])
param_df = pd.DataFrame([])
for g, d in grouped.groupby('gluc_conc_mM'):
    d.sort_values('time_hr', inplace=True)
    d = d[d['time_hr'] >= 2]
    _diff = np.round(np.diff(d['od_600nm_sub'].values), decimals=4)
    ind = np.where(_diff <= 0)[0][0]
    
    # Find the max time and the yield
    _time = d['time_hr'].values[ind-15]
    _od = d['od_600nm_sub'].values[ind - 1]

    # Assemble the yield df 
    yield_df = yield_df.append({'gluc_conc_mM':g, 'max_od_600nm':_od}, 
                                ignore_index=True)

    # crop the fitting data
    _d = d[d['time_hr'] <= _time] 
    crop_dfs.append(_d)

    # Compute the fit. 
    popt = scipy.stats.linregress(_d['time_hr'], np.log(_d['od_600nm_sub']))

    # Store the fit parameters
    param_df = param_df.append({'gluc_conc_mM':g, 
                                'growth_rate': popt[0],
                                'intercept': popt[1],
                                'growth_rate_err': popt[-1]},
                                ignore_index=True)
    time_range = np.linspace(2, _time, 100)
    _fit = np.exp(popt[1] + popt[0] * time_range)
    _fit_df = pd.DataFrame([])
    _fit_df['time_hr'] = time_range
    _fit_df['od_600nm_sub'] = _fit
    _fit_df['gluc_conc_mM'] = g
    fit_dfs.append(_fit_df)

    # Print the inferred growth rate
    print(f'Growth rate for {g} mM Glucose: λ = {popt[0]:0.3f} ± {popt[-1]:0.3f} per hr.')

crop_df = pd.concat(crop_dfs, sort=False)
fit_df = pd.concat(fit_dfs, sort=False)


points = alt.Chart(crop_df).mark_point(size=10, opacity=0.75).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density', scale=alt.Scale(type='log')),
        color=alt.Color('gluc_conc_mM:Q', title='glucose [mM]')
)

fits = alt.Chart(fit_df).mark_line(size=0.5).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density', scale=alt.Scale(type='log')),
        color=alt.Color('gluc_conc_mM:Q', title='glucose [mM]')
)

fits + points

# %%

# Given the y
# Plot the yield
yield_points = alt.Chart(yield_df).mark_point(size=80).encode(
        x=alt.X('gluc_conc_mM:Q', title='glucose [mM]'),
        y=alt.Y('max_od_600nm:Q', title='saturating OD [a.u.]')
)

yield_points
# %%
 