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
data = pd.read_csv('../../../../data/plate_reader/2021-04-19_NCM3722_glucose_tecan_bmg_comparison/processed/2021-04-19_NCM3722_glucose_tecan_bmg_comparison_tidy.csv')

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
)

save(raw, './output/2021-04-19_NCM3722_glucose_raw_traces.pdf')
raw

# Compute the average
# avg = subbed.groupby(['instrument', 'time_hr']).mean().reset_index()
# %%
base = alt.Chart(subbed).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
             color=alt.Color('instrument:N', title='instrument')
)

avg = base.mark_line().encode(
        y=alt.Y('mean(od_600nm_sub):Q', title='optical density [a.u.]')
)

band = base.mark_errorband(extent='stdev').encode(
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]')
       
)
band + avg
# %%

# Analyze only the spark data 
spark = subbed[subbed['instrument']=='Tecan Spark']
avg = spark.groupby('time_hr').mean().reset_index()

# Crop to the exponential region
crop = avg[avg['time_hr']<= 3.5]

popt  = scipy.stats.linregress(crop['time_hr'], np.log(crop['od_600nm_sub'].values))
print(f'Estimated λ ≈ {popt[0]:0.3f} ± {popt[-1]:0.3f} per hr')
time_range = np.linspace(0,4.1, 100)
fit = np.exp(popt[1] + popt[0] * time_range)
_df = pd.DataFrame([])
_df['time_hr'] = time_range
_df['od_600nm_sub'] = fit
points = alt.Chart(crop).mark_point(size=80, opacity=0.5).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]',
            scale=alt.Scale(zero=False)),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]',
                scale=alt.Scale(type='log'))
)

fit = alt.Chart(_df).mark_line().encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u]')
)

tecan_layer = (points + fit).properties(title=f'Tecan Spark; λ ≈ {popt[0]:0.3f} ± {popt[-1]:0.3f} per hr')

# %%

# Analyze only the spark data 
bmg = subbed[subbed['instrument']=='BMG Omega']
avg = bmg.groupby('time_hr').mean().reset_index()

# Crop to the exponential region
crop = avg[(avg['time_hr']<= 3.5) & (avg['time_hr'] >= 2)]

popt  = scipy.stats.linregress(crop['time_hr'], np.log(crop['od_600nm_sub'].values))
print(f'Estimated λ ≈ {popt[0]:0.3f} ± {popt[-1]:0.3f} per hr')

time_range = np.linspace(1.8,4, 100)
fit = np.exp(popt[1] + popt[0] * time_range)
_df = pd.DataFrame([])
_df['time_hr'] = time_range
_df['od_600nm_sub'] = fit
points = alt.Chart(crop).mark_point(size=80, opacity=0.5).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u.]',
                scale=alt.Scale(type='log'))
)

fit = alt.Chart(_df).mark_line().encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]',
            scale=alt.Scale(zero=False)),
        y=alt.Y('od_600nm_sub:Q', title='optical density [a.u]')
)

points + fit
bmg_layer = (points + fit).properties(title=f'BMG Omega; λ ≈ {popt[0]:0.1f} ± {popt[-1]:0.1f} per hr')

# %%
save(tecan_layer | bmg_layer, './output/2021-04-19_fit_comparison.pdf')

# %%

