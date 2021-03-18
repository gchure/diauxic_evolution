#%%
import numpy as np 
import pandas as pd 
import diaux.viz 
import scipy.stats
import altair as alt
from altair_saver import save
colors, palette = diaux.viz.altair_style()
alt.data_transformers.disable_max_rows()

# Load the tidied data 
data = pd.read_csv('../../../data/2021-03-15_REL606_NCM3722_growth_medium_comparison/processed/2021-03-15_REL606_NCM3722_NC_DM_comparison_tidy.csv')

#%%
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
for g, d in sub.groupby(['strain', 'medium', 'glucose_mM', 'replicate']):
    d.sort_values(by='time_s', inplace=True)
    d['od_600nm_rel'] = d['od_600nm_sub'].values / d['od_600nm_sub'].values[0]
    rel.append(d)
rel = pd.concat(rel)

# Convert time to hrs. 
rel['time_hr'] = rel['time_s'] / 3600

# Compute the mean and std_err for each. 
grouped = rel.groupby(['time_s', 'strain', 'glucose_mM', 'medium']).mean().reset_index()

# %%
# Set up the chart. 
base = alt.Chart(grouped)
points = base.mark_line().encode(
            x=alt.X(field='time_hr', type='quantitative', title = 'elapsed time [hr]'),
            y=alt.Y('od_600nm_sub:Q', 
                    title='optical density', scale=alt.Scale(type='log')),
            color=alt.Color('glucose_mM:O', title='glucose mM')
            ).facet(column='medium:N', row='strain:N')

save(points, './output/2021-03-15_REL606_NCM3722_growth_medium_comparison.pdf')
points

# %%
# Do a simple estimate for the growth rate between 1 and 2.5 hrs elapsed
early_exp = grouped[(grouped['time_hr'] >= 4) &
                    (grouped['time_hr'] <= 6)]

time_range = np.linspace(4, 6, 2020)
fit_df = []
for g, d in early_exp.groupby(['strain', 'medium', 'glucose_mM']):
    popt = scipy.stats.linregress(d['time_hr'], np.log(d['od_600nm_sub']))
    fit = np.exp(popt[0] * time_range + popt[1]) 
    _df = pd.DataFrame([])
    _df['od_600nm'] = fit
    _df['time'] = time_range
    _df['strain'] = g[0]
    _df['medium'] = g[1]
    _df['glucose_mM'] = g[2]
    fit_df.append(_df)
    print(f'{g[0]}, {g[1]}, {g[2]} mM glucose, λ = {popt[0]:0.3f} ± {popt[-1]:0.3f} / hr')
fit_df = pd.concat(fit_df, sort=False)

fit_base = alt.Chart(fit_df)
point_base = alt.Chart(early_exp)
fits = fit_base.mark_line().encode(
        x=alt.X('time:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm:Q', scale=alt.Scale(type='log'), title='optical density [a.u.]'),
        color=alt.Color('medium:N', title='growth medium')
        )

points= point_base.mark_point(size=10).encode(
        x=alt.X('time_hr:Q', title='elapsed time [hr]'),
        y=alt.Y('od_600nm_sub:Q', scale=alt.Scale(type='log'), title='optical density [a.u.]'),
        color=alt.Color('medium:N', title='growth medium'),
        )
# (points + fits).facet(row='strain:N', column='glucose_mM:N')
points

# %%
