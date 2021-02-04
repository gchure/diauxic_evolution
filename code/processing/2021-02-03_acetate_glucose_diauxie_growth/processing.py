#%%
import pandas as pd 
import numpy as np 
import scipy.stats
import altair as alt
from altair_saver import save
import diaux.viz 
colors, palette = diaux.viz.altair_style()

#%%
data = pd.read_csv('./input/2021-02-03_acetate_glucose_diauxie_growth.csv')

# Clean to the linear decade
data = data[(data['od_600nm']>=0.04) & (data['od_600nm'] <= 0.4)]

# Convert clock time to a datetime  
data['clock_time'] = pd.to_datetime(data['clock_time'])

# # Compute the elapsed time
dfs = []
for g, d in data.groupby(['media']):
    d.sort_values(by='clock_time', inplace=True) 
    d['elapsed_time_hr'] = d['clock_time'] - d['clock_time'].values[0]
    d['rel_od_600nm'] = d['od_600nm'].values / d['od_600nm'].values[0]
    dfs.append(d)   
data = pd.concat(dfs, sort=False)
data['elapsed_time_hr'] = (data['elapsed_time_hr'].astype('timedelta64[m]'))/60


# For the specified media, compute the growth rate. 
time = np.linspace(0, 5, 200)
fit_dfs = []
for g, d in data.groupby(['media']):
    slope, inter, _, _, se = scipy.stats.linregress(d['elapsed_time_hr'].values,
                                  np.log(d['od_600nm'].values))
    fit = np.exp(inter + slope * time )  
    label  = f'λ = {slope:0.03f} ± {se:0.03f} / hr'
    print(f'{g} medium; {label}')
    fit_df = pd.DataFrame([])
    fit_df['elapsed_time_hr'] = time
    fit_df['od_600nm'] = fit
    fit_df['rel_od_600nm'] = fit / fit[0]
    fit_df['media'] = g
    fit_df['label'] = label
    fit_dfs.append(fit_df)
fit_df = pd.concat(fit_dfs, sort=False)


# %%
# Set up the chart base
od_points = alt.Chart(data).mark_point(size=80).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='od_600nm', type='quantitative', title='optical density [a.u.]',
                    scale=alt.Scale(type='log', domain=[0.04, 0.4])),
            color=alt.Color(field='media', type='nominal', title='medium type'))

od_fit = alt.Chart(fit_df).mark_line(clip=True).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='od_600nm', type='quantitative', title='optical density [a.u.]',
                    scale=alt.Scale(type='log', domain=[0.04, 0.4])),
            color=alt.Color(field='label', type='nominal', title='growth rate',
                            scale=alt.Scale(domain=fit_df['label'].unique(), range=palette[:4],)))

od = (od_fit + od_points).resolve_scale(color='independent')

rel_od_points = alt.Chart(data).mark_point(size=80).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='rel_od_600nm', type='quantitative', title='relative optical density',
                    scale=alt.Scale(type='log', domain=[1, 10])),
            color=alt.Color(field='media', type='nominal', title='medium type'))

rel_od_fit = alt.Chart(fit_df).mark_line(clip=True).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='rel_od_600nm', type='quantitative', title='relative optical density',
                    scale=alt.Scale(type='log', domain=[1, 10])),
            color=alt.Color(field='label', type='nominal', title='growth rate',
                            scale=alt.Scale(domain=fit_df['label'].unique(), range=palette[:4],)))

layer = (od_fit + od_points).resolve_scale(color='independent') & (rel_od_fit + rel_od_points).resolve_scale(color='independent')
layer
# save(layer, './output/2021-02-03_acetate_glucose_diauxie.pdf')

# %%
