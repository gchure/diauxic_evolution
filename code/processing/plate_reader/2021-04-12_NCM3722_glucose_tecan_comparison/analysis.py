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
DATA_PATH = '../../../../data/plate_reader/2021-04-12_NCM3722_glucose_tecan_comparison/processed/'
data = pd.read_csv(f'{DATA_PATH}/2021-04-12_NCM3722_tecan_comparison_tidy.csv')

# Convert time to hours 
data['time_hr'] = data['time_s'] / 3600
# Compute the replicate averages
# grouped = data.groupby(['instrument', 'strain', 'time_hr', 'time_s'])['od_600nm'].agg(('mean', 'std')).reset_index()

# Compute the sub
dfs = []
for g, d in data.groupby(['instrument', 'time_s', 'time_hr']):
    avg_blank = d[d['strain']=='blank']['od_600nm'].mean()
    d['od_600nm_sub'] = d['od_600nm'].values -  avg_blank
    dfs.append(d)
subbed = pd.concat(dfs, sort=False)
subbed = subbed[subbed['strain']!='blank']
subbed
# %%
raw = alt.Chart(subbed).mark_line().encode(
    x=alt.X('time_hr:Q', title='elapsed time [hr]'),
    y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]',
            scale=alt.Scale(type='log')),
    color=alt.Color('instrument:N', title='tecan instrument'),
    opacity=alt.Opacity('replicate:O', legend=None)
)

save(raw, 'output/2021-04-12_tecan_comparison_raw.pdf')
# %%

# Plot the averages
avg = alt.Chart(subbed).mark_line().encode(
    x=alt.X('time_hr:Q', title='elapsed time [hr]'),
    y=alt.Y('mean(od_600nm_sub):Q', title='optical density [a.u.]'),
    color=alt.Color('instrument:N', title='tecan instrument'))



err = alt.Chart(subbed).mark_errorband(extent='stdev').encode(
    x=alt.X('time_hr:Q', title='elapsed time [hr]'),
    y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]'),
    color=alt.Color('instrument:N', title='tecan instrument'))


layer = avg + err

save(layer, 'output/2021-04-12_tecan_comparison_avg.pdf')
layer


# %%
# Consider only the first five hours
crop = subbed[subbed['time_hr'] <= 5]

# Compute the average
crop = crop.groupby(['instrument', 'time_hr']).mean().reset_index()

# for each instrument, to the regression
params = {}
fits = []
time_range = np.linspace(0, 5, 200)
for g, d in crop.groupby(['instrument']):
    popt = scipy.stats.linregress(d['time_hr'], np.log(d['od_600nm_sub']))
    params[g] = {'slope': popt[0], 'inter':popt[1], 'err':popt[-1]}
    fit = np.exp(popt[1] + popt[0] * time_range)
    _df = pd.DataFrame([])
    _df['time_hr'] = time_range 
    _df['od_600nm_sub'] = fit
    _df['instrument'] = g
    fits.append(_df)
    print(f'Estimate for Tecan {g}: λ ≈ {popt[0]:0.3f} ± {popt[-1]:0.3f} per hr')
fits = pd.concat(fits, sort=False)

# %%
points = alt.Chart(crop).mark_point(size=50, opacity=0.5).encode(
    x=alt.X('time_hr:Q', title='elapsed time [hr]'),
    y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]',
            scale=alt.Scale(type='log')),
    color=alt.Color('instrument:N', title='tecan instrument')
)
fit = alt.Chart(fits).mark_line().encode(
    x=alt.X('time_hr:Q', title='elapsed time [hr]'),
    y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]',
            scale=alt.Scale(type='log')),
    color=alt.Color('instrument:N', title='tecan instrument')
)

layer = (points + fit).resolve_scale(color='independent')
save(layer, './output/2021-04-12_tecan_comparison_fits.pdf')
# %%
