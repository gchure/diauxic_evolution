#%%
import pandas as pd 
import numpy as np 
import scipy.stats
import altair as alt
from altair_saver import save
import diaux.viz 
colors, palette = diaux.viz.altair_style()

data = pd.read_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/raw/2021-02-15_NCM3722_growth_curve.csv')
hplc = pd.read_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_diauxic_shift_glucose_acetate_turnover.csv')

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

# %% Merge the HPLC data on the growth curve
for g, d in data.groupby(['sample_id']):
    # get the correct hplc point
    if g != 0:
        hplc_point = hplc[hplc['sample_id']==g]

        # Get the glucose and acetate concentrations
        glucose = hplc_point[hplc_point['carbon_source']=='glucose']['rescaled_concentration_mM'].values[0]
        acetate = hplc_point[hplc_point['carbon_source']=='acetate']['rescaled_concentration_mM'].values[0]

        # locate and insert the concentration information
        data.loc[data['sample_id']==g, 'glucose_concentration_mM'] = glucose
        data.loc[data['sample_id']==g, 'acetate_concentration_mM'] = acetate
#%% Save the data to disk
data.to_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_glucose_acetate_turnover_growth.csv', index=False)

#%%
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

#%%
rel_od_points = alt.Chart(data).mark_point(size=90, opacity=0.75).encode(
    x=alt.X(field='elapsed_time_hr',type='quantitative', title='elapsed time [hr]',
            scale=alt.Scale(domain=[0, 7])),
    y=alt.Y(field='rel_od_600nm', type='quantitative', title='optical density [a.u.]',
            scale=alt.Scale(type='log')),
    color=alt.Color(field='region', type='nominal', title='growth phase'),
    shape=alt.Shape(field='sampled', type='nominal', title='HPLC sampled?')
    )

rel_od_points


save(rel_od_points, f'output/2021-02-15_diauxic_shift_curve.pdf')
#%%


# %%
