#%%
import numpy as np 
import pandas as pd 
import diaux.viz 
import scipy.stats
import altair as alt
from altair_saver import save
colors, palette = diaux.viz.altair_style()

# Load the tidied data 
data = pd.read_csv('./output/2021-02-11_NCM_REL_diauxie_tidy.csv')

# Do the background subtraction
sub = []
for g, d in data.groupby(['time_s']):
    # Compute the average of all of the blanks
    blank = d[d['strain']=='blank']['od_600nm'].mean()
    d['od_600nm_sub'] = d['od_600nm'].values - blank
    sub.append(d)
sub = pd.concat(sub, sort=False)
sub = sub[sub['strain'] != 'blank']
# Compute the relative change in OD
rel = []
for g, d in sub.groupby(['strain', 'medium', 'replicate']):
    d.sort_values(by='time_s', inplace=True)
    d['od_600nm_rel'] = d['od_600nm_sub'].values / d['od_600nm_sub'].values[0]
    rel.append(d)
rel = pd.concat(rel)

# Convert time to hrs. 
rel['time_hr'] = rel['time_s'] / 3600

# Compute the mean and std_err for each. 
grouped = rel.groupby(['time_s', 'strain', 'medium']).mean().reset_index()

# %%
# Set up the chart. 
base = alt.Chart(grouped)
points = base.mark_line().encode(
            x=alt.X(field='time_hr', type='quantitative', title = 'elapsed time [hr]'),
            y=alt.Y(field='od_600nm_sub', type='quantitative', 
                    title='optical density', scale=alt.Scale(type='log')),
            color=alt.Color(field='strain', type='nominal', title='strain')
            ).facet(column='medium:N')

save(points, './output/2021-02-11_NCM_REL_diauxie_lagtime_absolute.pdf')
points

# %%
# Do a simple estimate for the growth rate between 1 and 2.5 hrs elapsed
early_exp = grouped[(grouped['time_hr'] >= 1) &
                    (grouped['time_hr'] <= 2.5)]

time_range = np.linspace(1, 2.5, 2020)
fit_df = []
for g, d in early_exp.groupby(['strain', 'medium']):
    popt = scipy.stats.linregress(d['time_hr'], np.log(d['od_600nm_sub']))
    print(f'{g[0]}, {g[1]}, λ = {popt[0]:0.3f} ± {popt[-1]:0.3f} / hr')

# %%
