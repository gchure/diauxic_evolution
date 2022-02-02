#%%
import numpy as np 
import pandas as pd 
import scipy.stats
import diaux.viz 
import altair as alt 
from altair_saver import save
alt.data_transformers.disable_max_rows()
colors, paleette = diaux.viz.altair_style()

# Load the tidy dataset
data = pd.read_csv('../../../../data/plate_reader/2021-04-27_NCM3722_glucose_acetate_tecan_biotek_comparison/processed/2021-04-27_NCM3722_glucose_acetate_tecan_biotek_comparison.csv')

# Convert time to hours 
data['time_hr'] = data['time_s'] / 3600
data.drop(columns='time_s', inplace=True)
# Compute the sub
dfs = []
for g, d in data.groupby(['instrument', 'time_hr']):
    avg_blank = d[d['sample']=='blank']['od_600nm'].mean()
    d['od_600nm_sub'] = d['od_600nm'].values -  avg_blank
    dfs.append(d)
subbed = pd.concat(dfs, sort=False)
subbed = subbed[subbed['sample']!='blank']
# %%
raw = alt.Chart(subbed).mark_line(size=1).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]'),
        color=alt.Color('instrument:N', title='instrument'),
        opacity=alt.Size('replicate:O', title='replicate', legend=None)
).facet(column='medium')

save(raw, './output/2021-04-27_NCM3722_glucose_acetate_raw_traces.pdf')
raw

# %%

# Generate a plot the avgs
base = alt.Chart(subbed).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
             color=alt.Color('instrument:N', title='instrument')
)

avg = base.mark_line().encode(
        y=alt.Y('mean(od_600nm_sub):Q', title='optical density [a.u.]',
            scale=alt.Scale(type='log'))
)

band = base.mark_errorband(extent='stdev').encode(
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]')
       
)
averaged = (band + avg).facet(column='medium')
save(averaged, './output/2021-04-27_NCM3722_glucose_acetate_avg_traces.pdf')
averaged
# %%

# Truncate and compute the inferred growth rates
gluc_bounds = [0, 2.6]
ac_bounds = [0, 6]
fit_dfs = []
crop_dfs = []
# grouped = subbed.groupby(['time_hr', 'medium', 'instrument']).mean().reset_index()
grouped = subbed
for medium, bounds in zip(['glucose', 'acetate'], [gluc_bounds, ac_bounds]):
    cropped = grouped[(grouped['time_hr']>=bounds[0]) & 
                      (grouped['time_hr'] <= bounds[1]) &
                      (grouped['medium']==medium)]
    crop_dfs.append(cropped)
    for g, d in cropped.groupby(['instrument']): 
        popt = scipy.stats.linregress(d['time_hr'], np.log(d['od_600nm_sub']))
        slope = popt[0]
        inter = popt[1]
        err = popt[-1]
        print(f'{g} growth in {medium}: λ = {slope:0.3} ± {err:0.3f} per hr')
        time_range = np.linspace(bounds[0], bounds[1], 200)
        best_fit = np.exp(inter + slope * time_range)
        _df = pd.DataFrame([])
        _df['time_hr'] = time_range 
        _df['od_600nm_sub'] = best_fit
        _df['medium'] = medium
        _df['instrument'] = g
        fit_dfs.append(_df)

fit_df = pd.concat(fit_dfs, sort=False)
crop_df = pd.concat(crop_dfs, sort=False)
# %%
# Set up the fit plots 
out = []
for m in ['glucose', 'acetate']:
    _crop = crop_df[crop_df['medium']==m]
    _fit = fit_df[fit_df['medium']==m]
    points = alt.Chart(_crop).mark_point(size=50, opacity=0.75).encode(
            x=alt.X('time_hr:Q', title='elapsed time [hr]'),
            y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]', 
                    scale=alt.Scale(type='log')),
            color=alt.Color('instrument:N', title='instrument')
            )

    fit = alt.Chart(_fit).mark_line().encode(
            x=alt.X('time_hr:Q', title='elapsed time [hr]'),
            y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]', 
                    scale=alt.Scale(type='log')),
            color=alt.Color('instrument:N', title='instrument')
            )
    _out = (points + fit).properties(title=f'{m} medium')
    out.append(_out)

out = out[0] | out[1] 
save(out, './output/2021-04-27_glucose_acetate_fit.pdf')
# %%
