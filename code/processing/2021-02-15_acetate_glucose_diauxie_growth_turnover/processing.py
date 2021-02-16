#%%
import pandas as pd 
import numpy as np 
import scipy.stats
import altair as alt
from altair_saver import save
import diaux.viz 
colors, palette = diaux.viz.altair_style()

data = pd.read_csv('./input/2021-02-15_NCM_acetate_glucose_diauxie_growth.csv')

# Clean to the linear decade
data = data[(data['od_600nm']>=0.04) & (data['od_600nm'] <= 0.405)]

# Add identifier of what is bing sampled
data['sampled'] = False
data.loc[data['sample_id'] > 0, 'sampled'] = True

# Convert clock time to a datetime  
data['clock_time'] = pd.to_datetime(data['clock_time'])

# Compute the elapsed time
data.sort_values(by='clock_time', inplace=True) 
data['elapsed_time_hr'] = data['clock_time'] - data['clock_time'].values[0]
data['rel_od_600nm'] = data['od_600nm'].values / data['od_600nm'].values[0]
data['elapsed_time_hr'] = (data['elapsed_time_hr'].astype('timedelta64[m]'))/60

# Add regional identifier -- this is roughly chosen by inspection
data.loc[data['elapsed_time_hr'] < 2.1, 'region'] = 'growth on glucose'
data.loc[(data['elapsed_time_hr'] >= 2.1) & (data['elapsed_time_hr'] < 5), 
        'region'] = 'lag phase'
data.loc[data['elapsed_time_hr'] >= 5, 'region'] = 'growth on acetate'

# For the specified media, compute the growth rate. 
time = np.linspace(0, 5, 200)
glucose_data = data[data['elapsed_time_hr'] <= 2]

fit_dfs = []
for g, d in data.groupby(['region']):
    slope, inter, _, _, se = scipy.stats.linregress(d['elapsed_time_hr'].values,
                                      np.log(d['od_600nm'].values))
    fit = np.exp(inter + slope * time )  
    label  = f'λ = {slope:0.03f} ± {se:0.03f} / hr'
    print(f'{g},{label}')
    fit_df = pd.DataFrame([])
    fit_df['elapsed_time_hr'] = time
    fit_df['od_600nm'] = fit
    fit_df['rel_od_600nm'] = fit / fit[0]
    fit_df['label'] = label
    fit_df['region'] = g   
    fit_dfs.append(fit_df)
fit_df = pd.concat(fit_dfs, sort=False)

rel_od_points = alt.Chart(data).mark_point(size=90, opacity=0.75).encode(
    x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]',
            scale=alt.Scale(domain=[0, 7])),
    y=alt.Y(field='rel_od_600nm', type='quantitative', title='optical density [a.u.]',
            scale=alt.Scale(type='log')),
    color=alt.Color(field='sampled', type='nominal', title='sampled collected'),
    # shape=alt.Shape(field='sampled', type='nominal', title='HPLC sampled?')
    )

# rel_od_fit = alt.Chart(fit_df).mark_line(clip=True, color=palette[0]).encode(
#     x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]'),
#     y=alt.Y(field='od_600nm', type='quantitative', title='optical density [a.u.]',
#             scale=alt.Scale(type='log')))


rel_od_points
# layer = (rel_od_points + rel_od_fit).resolve_scale(color='independent')

save(rel_od_points, f'output/relative_od_shift.pdf')

# %%
