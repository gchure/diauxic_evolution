#%%
import pandas as pd 
import numpy as np 
import scipy.stats
import altair as alt
from altair_saver import save
import diaux.viz 
colors, palette = diaux.viz.altair_style()

#%%
data = pd.read_csv('./input/2021-02-08_NCM_REL_diauxie_growth.csv')

# Clean to the linear decade
data = data[(data['od_600nm']>=0.04) & (data['od_600nm'] <= 0.4)]

# Convert clock time to a datetime  
data['clock_time'] = pd.to_datetime(data['clock_time'])

# # Compute the elapsed time
dfs = []
for g, d in data.groupby(['media', 'strain']):
    d.sort_values(by='clock_time', inplace=True) 
    d['elapsed_time_hr'] = d['clock_time'] - d['clock_time'].values[0]
    d['rel_od_600nm'] = d['od_600nm'].values / d['od_600nm'].values[0]
    dfs.append(d)   
data = pd.concat(dfs, sort=False)
data['elapsed_time_hr'] = (data['elapsed_time_hr'].astype('timedelta64[m]'))/60


# For the specified media, compute the growth rate. 
time = np.linspace(0, 5, 200)
fit_dfs = []
for g, d in data.groupby(['media', 'strain']):
    slope, inter, _, _, se = scipy.stats.linregress(d['elapsed_time_hr'].values,
                                  np.log(d['od_600nm'].values))
    fit = np.exp(inter + slope * time )  
    label  = f'λ = {slope:0.03f} ± {se:0.03f} / hr'
    if g[0] == 'shift':
            continue
    print(f'strain {g[1]}, {g[0]}; {label}')
    fit_df = pd.DataFrame([])
    fit_df['elapsed_time_hr'] = time
    fit_df['od_600nm'] = fit
    fit_df['rel_od_600nm'] = fit / fit[0]
    fit_df['media'] = g[0]
    fit_df['strain'] = g[1]
    fit_df['label'] = label
    fit_dfs.append(fit_df)
fit_df = pd.concat(fit_dfs, sort=False)

# Make different plots for different strains
for g, d in data.groupby(['media']):
        fit = fit_df[fit_df['media']==g]
        rel_od_points = alt.Chart(d).mark_point(line=True, size=80, opacity=0.75).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='rel_od_600nm', type='quantitative', title='relative optical density',
                    scale=alt.Scale(type='log', domain=[1, 10])),
            color=alt.Color(field='strain', type='nominal', title='strain')
            ).properties(title=f'growth in {g}')

        rel_od_fit = alt.Chart(fit).mark_line(clip=True).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='rel_od_600nm', type='quantitative', title='relative optical density',
                    scale=alt.Scale(type='log', domain=[1, 10])),
            color=alt.Color(field='strain', type='nominal'))

        if len(fit) > 0:
                layer = (rel_od_points + rel_od_fit).resolve_scale(color='independent')
        else:
                layer = rel_od_points
        save(layer, f'output/{g}_medium_relative_od.pdf')

# %%
# Set up the chart base
od_points = alt.Chart(data).mark_point(size=80).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='od_600nm', type='quantitative', title='optical density [a.u.]',
                    scale=alt.Scale(type='log', domain=[0.04, 0.4])),
            shape=alt.Shape(field='strain', type='nominal', title='strain'),
            color=alt.Color(field='media', type='nominal', title='medium type'))

od_fit = alt.Chart(fit_df).mark_line(clip=True).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='od_600nm', type='quantitative', title='optical density [a.u.]',
                    scale=alt.Scale(type='log', domain=[0.04, 0.4])),   
            strokeDash=alt.StrokeDash(field='strain', type='nominal', title='strain'),
            color=alt.Color(field='label', type='nominal', title='growth rate',
                            scale=alt.Scale(domain=fit_df['label'].unique(), range=palette[:4],)))

od = (od_fit + od_points).resolve_scale(color='independent')

rel_od_points = alt.Chart(data).mark_point(size=80, opacity=0.75).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='rel_od_600nm', type='quantitative', title='relative optical density',
                    scale=alt.Scale(type='log', domain=[1, 10])),
            shape=alt.Shape(field='strain', type='nominal', title='strain'),
            color=alt.Color(field='media', type='nominal', title='medium type'))

rel_od_fit = alt.Chart(fit_df).mark_line(clip=True).encode(
            x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
            y=alt.Y(field='rel_od_600nm', type='quantitative', title='relative optical density',
                    scale=alt.Scale(type='log', domain=[1, 10])),
            strokeDash=alt.StrokeDash(field='strain', type='nominal', legend=None),
            color=alt.Color(field='media', type='nominal', legend=None))
                            

layer = (od_fit + od_points).resolve_scale(color='independent') & (rel_od_fit + rel_od_points).resolve_scale(color='independent')
layer
# save(layer, './output/2021-02-03_acetate_glucose_diauxie.pdf')

        # %%
