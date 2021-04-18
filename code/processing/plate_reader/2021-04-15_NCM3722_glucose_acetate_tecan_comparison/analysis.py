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
DATA_PATH = '../../../../data/plate_reader/2021-04-15_NCM3722_glucose_acetate_tecan_comparison/processed/'
data = pd.read_csv(f'{DATA_PATH}/2021-04-15_NCM3722_glucose_acetate_tecan_comparison_tidy.csv')

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
subbed
# %%
raw = alt.Chart(subbed).mark_line().encode(
    x=alt.X('time_hr:Q', title='elapsed time [hr]'),
    y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]',
            scale=alt.Scale(type='log')),
    color=alt.Color('instrument:N', title='tecan instrument'),
    opacity=alt.Opacity('replicate:O', legend=None)
).facet(column='media')
raw
save(raw, 'output/2021-04-15_glucose_acetate_tecan_comparison_raw.pdf')
# %%

# Plot the averages
out = []
for media in ['glucose', 'acetate']:
    _data = subbed[subbed['media']==media]
    avg = alt.Chart(_data).mark_line().encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('mean(od_600nm_sub):Q', title='optical density [a.u.]'),
        # scale=alt.Scale(type='log')),
        color=alt.Color('instrument:N', title='tecan instrument'))

    err = alt.Chart(_data).mark_errorband(extent='stdev').encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]'),
        color=alt.Color('instrument:N', title='tecan instrument'))


    layer = (avg + err).properties(title=f'N-C- + {media}')
    out.append(layer)

result = out[0] | out[1]
result
save(result, 'output/2021-04-15_glucose_acetate_tecan_comparison_avg.pdf')



# %%
# Consider only the first five hogurs
gluc_bounds = [1, 4]
ac_bounds = [3, 6]


fits = []
for media, bounds in zip(['glucose', 'acetate'], [gluc_bounds, ac_bounds]):
    crop = subbed[(subbed['media']==media) & ((subbed['time_hr'] >= bounds[0]) &
                   (subbed['time_hr'] <= bounds[1]))]

    # Compute the average
    crop = crop.groupby(['instrument', 'time_hr']).mean().reset_index()

    # for each instrument, to the regression
    params = {}
    time_range = np.linspace(bounds[0] - 0.5, bounds[1] + 0.5, 200)
    for g, d in crop.groupby(['instrument']):
        popt = scipy.stats.linregress(d['time_hr'], np.log(d['od_600nm_sub']))
        params[g] = {'media': media, 'slope': popt[0], 'inter':popt[1], 'err':popt[-1]}
        fit = np.exp(popt[1] + popt[0] * time_range)
        _df = pd.DataFrame([])
        _df['time_hr'] = time_range 
        _df['od_600nm_sub'] = fit
        _df['instrument'] = g
        _df['media'] = media
        fits.append(_df)
        print(f'Estimate for Tecan {g}, {media} medium: λ ≈ {popt[0]:0.3f} ± {popt[-1]:0.3f} per hr')        

fits = pd.concat(fits, sort=False)


# %%
out = []
for media, bounds in zip(['glucose', 'acetate'], [gluc_bounds, ac_bounds]):
    crop = subbed[(subbed['media']==media) & ((subbed['time_hr'] >= bounds[0]) &
                   (subbed['time_hr'] <= bounds[1]))]

    # Compute the average
    crop = crop.groupby(['instrument', 'time_hr']).mean().reset_index()
    _fits = fits[fits['media']==media]

    points = alt.Chart(crop).mark_point(size=50, opacity=0.5).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]',
            scale=alt.Scale(type='log')),
        color=alt.Color('instrument:N', title='tecan instrument')
    )

    fit = alt.Chart(_fits).mark_line().encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]',
            scale=alt.Scale(type='log')),
        color=alt.Color('instrument:N', title='tecan instrument')
    )

    layer = (points + fit).resolve_scale(color='independent').properties(title='N-C- + {media}')
    out.append(layer)
layer = out[0] | out[1]
layer
save(layer, './output/2021-04-15_glucose_acetate_tecan_comparison_fits.pdf')
# %%
out = []
for media in ['glucose', 'acetate']:
    _data = subbed[subbed['media']==media]
    avg = alt.Chart(_data).mark_line(clip=True).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]',
                scale=alt.Scale(domain=[0, 6])),
        y=alt.Y('mean(od_600nm_sub):Q', title='optical density [a.u.]'),
        # scale=alt.Scale(type='log')),
        color=alt.Color('instrument:N', title='tecan instrument'))

    err = alt.Chart(_data).mark_errorband(extent='stdev', clip=True).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]'),
        color=alt.Color('instrument:N', title='tecan instrument'))


    layer = (avg + err).properties(title=f'N-C- + {media}')
    out.append(layer)

result = out[0] | out[1]
result
# save(result, 'output/2021-04-15_glucose_acetate_tecan_comparison_avg_crop.pdf')



# %%
