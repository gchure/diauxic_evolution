#%%
sys.path.insert(0, '../../../../')
import numpy as np 
import pandas as pd 
import altair as alt 
from altair_saver import save
import scipy.stats
import diaux.viz 
colors, palette = diaux.viz.altair_style()


# %%
# Load the various data sets
DATA_PATH = '../../../../data/metabolite_turnover/2021-04-04_REL606_glucose_turnover/processed/'
peaks = pd.read_csv(f'{DATA_PATH}/2021-04-04_REL606_glucose_peak_table.csv')
cal_peaks = pd.read_csv('../../../../data/hplc_calibration/2021-04-05_NC_DM_calibration/processed/2021-04-05_NC_DM_carbon_source_calibration.csv')
cal_peaks = cal_peaks[cal_peaks['buffer_base']=='DM']
growth = pd.read_csv(f'{DATA_PATH}/2021-04-04_REL606_glucose_growth.csv')

#%%
# For the peaks, map time point to OD
dfs = []
for g, d in peaks.groupby(['replicate', 'time_idx']):
    _growth = growth[(growth['replicate']==g[0]) & 
                     (growth['time_idx']==g[1])]['od_600nm'].values[0]
    d['od_600nm'] = _growth
    dfs.append(d)
peaks = pd.concat(dfs, sort=False)

# %%
# Make a plot of all of the compound areas as a function of time point
chart = alt.Chart(peaks, width=250, height=250
                 ).mark_point(size=80, opacity=0.75).encode(
            x=alt.X('time_idx:O', title='time point'),
            y=alt.Y('area:Q', title='integrated signal [mV]'),
            color=alt.Color('compound:N', title='compound')
).facet(column='replicate:N')
save(chart, './output/2021-04-04_REL606_compound_area_variation.png')
# %%
# Normalize the peak areas to phosphate
rel_peaks = []
for g, d in peaks.groupby(['replicate', 'time_idx']):
    phos_peak = d[d['compound']=='phosphate']['area'].values[0]
    d['rel_area_phosphate'] = d['area'].values / phos_peak
    rel_peaks.append(d)
rel_peaks = pd.concat(rel_peaks, sort=False)
# %%
chart = alt.Chart(rel_peaks, width=250, height=250
                 ).mark_point(size=80, opacity=0.75).encode(
            x=alt.X('time_idx:O', title='time point'),
            y=alt.Y('rel_area_phosphate:Q', 
                    title='signal relative to phosphate'),
            color=alt.Color('compound:N', title='compound')
).facet(column='replicate:N')
save(chart, './output/2021-04-04_REL606_compound_area_variation_normalized.png')

#%% 
# Set up the calibration curve for the data relative to the phosphate peak
rel_cal_peaks = []
for g, d in cal_peaks.groupby(['carbon_conc_mM']):
    phos_peak = d[d['compound']=='phosphate']['area'].values[0]
    d['rel_area_phosphate'] = d['area'].values / phos_peak
    rel_cal_peaks.append(d)
rel_cal_peaks = pd.concat(rel_cal_peaks, sort=False)
rel_cal_peaks = rel_cal_peaks[rel_cal_peaks['compound'].isin(['glucose', 'acetate'])]


# Perform the estimate
conc_range = np.linspace(0, 30, 100)
params = {}
fit_dfs = []
for g, d in rel_cal_peaks.groupby(['compound']):
    popt = scipy.stats.linregress(d['carbon_conc_mM'], d['rel_area_phosphate'])
    slope = popt[0]
    intercept = popt[1]
    params[g] = {'slope':slope, 'intercept':intercept}
    fit_df = pd.DataFrame([])
    fit_df['carbon_conc_mM'] = conc_range
    fit_df['rel_area_phosphate'] = intercept + slope * conc_range
    fit_df['compound'] = g
    fit_dfs.append(fit_df)
fit_df = pd.concat(fit_dfs, sort=False)


points = alt.Chart(rel_cal_peaks).mark_point(size=80).encode(
            x='carbon_conc_mM:Q',
            y='rel_area_phosphate:Q',
            color='compound:N'
)
fit = alt.Chart(fit_df).mark_line(size=2).encode(
            x=alt.X('carbon_conc_mM:Q', title='concentration [mM]'),
            y=alt.Y('rel_area_phosphate:Q', title='signal relative to phosphate'),
            color=alt.Color('compound:N')
)

layer = fit + points
save(layer, './output/DM_calibration_curve.png')
# %%
# Look at only glucose and acetate
samp_data = rel_peaks[rel_peaks['compound'].isin(['glucose', 'acetate'])]

glucose_conc = (samp_data[samp_data['compound']=='glucose']['rel_area_phosphate'] -\
                params['glucose']['intercept']) / params['glucose']['slope']
acetate_conc = (samp_data[samp_data['compound']=='acetate']['rel_area_phosphate'] -\
                params['acetate']['intercept']) / params['acetate']['slope']

# Conver the relative area to concentration measurements
samp_data.loc[samp_data['compound']=='glucose', 'conc_mM'] = glucose_conc
samp_data.loc[samp_data['compound']=='acetate', 'conc_mM'] = acetate_conc


# Do a simple linear regression of the yields, pooling all replicates 
yield_params = {}
od_range = np.linspace(0.005, 0.45, 200)
fit_df = []
for g, d in samp_data.groupby(['compound']):
    popt = scipy.stats.linregress(d['od_600nm'], d['conc_mM'])
    slope, intercept = popt[:2]
    err = popt[-1]
    yield_params[g] = {'slope':slope, 'inter':intercept, 'err':err}
    _fit_df = pd.DataFrame([])
    _fit_df['od_600nm'] = od_range
    _fit_df['conc_mM_lower'] = intercept + (slope-err) * od_range
    _fit_df['conc_mM_upper'] = intercept + (slope+err) * od_range
    _fit_df['conc_mM'] = intercept + (slope) * od_range
    _fit_df['compound'] = g
    fit_df.append(_fit_df)
fit_df = pd.concat(fit_df, sort=False)

#%%
ac_data = samp_data[samp_data['compound']=='acetate']
ac_fit = fit_df[fit_df['compound']=='acetate']
glu_data = samp_data[samp_data['compound']=='glucose']
glu_fit = fit_df[fit_df['compound']=='glucose']

#%%
glu_points = alt.Chart(glu_data, width=350, height=300
                      ).mark_point(size=80, opacity=0.95
                      ).encode(
                          x=alt.X('od_600nm:Q', title='optical density [a.u.]'),
                          y=alt.Y('conc_mM:Q', title='glucose concentration [mM]',
                                  scale=alt.Scale(domain=[5, 12])),
                          color=alt.Color('replicate:N', title='biological replicate')
                      )
glu_fit_base = alt.Chart(glu_fit, width=350, height=300
                      ).encode(
                          x=alt.X('od_600nm:Q', title='optical density [a.u.]'))
                                              
glu_uncertainty = glu_fit_base.mark_area(opacity=0.25).encode(
                        y=alt.Y('conc_mM_lower:Q', title='glucose concentration [mM]',
                                  scale=alt.Scale(domain=[7, 11])),
                          y2='conc_mM_upper:Q')
glu_best = glu_fit_base.mark_line(size=2).encode(
                        y=alt.Y('conc_mM:Q', title='glucose concentration [mM]',
                                  scale=alt.Scale(domain=[7, 11])))


glu_layer = (glu_uncertainty + glu_best + glu_points).properties(
                title=f"glucose yield = {yield_params['glucose']['slope']:0.1f} ± {yield_params['glucose']['err']:0.1f} mM / OD")
glu_layer

ac_points = alt.Chart(ac_data, width=350, height=300
                      ).mark_point(size=80, opacity=0.95
                      ).encode(
                          x=alt.X('od_600nm:Q', title='optical density [a.u.]'),
                          y=alt.Y('conc_mM:Q', title='acetate concentration [mM]',
                                  scale=alt.Scale(domain=[0, 1.5])),
                          color=alt.Color('replicate:N', title='biological replicate')
                      )
ac_fit_base = alt.Chart(ac_fit, width=350, height=300
                      ).encode(
                          x=alt.X('od_600nm:Q', title='optical density [a.u.]'))
                                              
ac_uncertainty = ac_fit_base.mark_area(opacity=0.25).encode(
                        y=alt.Y('conc_mM_lower:Q', title='acetate concentration [mM]',
                                  scale=alt.Scale(domain=[0, 1.5])),
                          y2='conc_mM_upper:Q')
ac_best = ac_fit_base.mark_line(size=2).encode(
                        y=alt.Y('conc_mM:Q', title='acetate concentration [mM]',
                                  scale=alt.Scale(domain=[0, 1.5])))


ac_layer = (ac_uncertainty + ac_best + ac_points).properties(
                title=f"acetate secretion = {yield_params['acetate']['slope']:0.2f} ± {yield_params['acetate']['err']:0.2f} mM / OD")

save(glu_layer | ac_layer, './output/2021-04-04_REL606_glucose_turnover.pdf')
save(glu_layer | ac_layer, './output/2021-04-04_REL606_glucose_turnover.png')

#%%
points = alt.Chart(samp_data, width=350, height=300).mark_point(size=80).encode(
            x=alt.X('od_600nm:Q', title='optical density [a.u.]'),
            y=alt.Y('conc_mM:Q', title='concentration [mM]'),
            color=alt.Color('replicate:N', title='biological replicate'),
            facet=alt.Facet('compound:N', header=alt.Header(labelFontSize=15)))

fit = alt.Chart(fit_df, width=350, height=300).mark_line(color=colors['black']).encode(
            x=alt.X('od_600nm:Q', title='optical density [a.u.]'),
            y=alt.Y('conc_mM:Q', title='concentration [mM]'),
            facet=alt.Facet('compound:N', header=alt.Header(labelFontSize=15)))

points + fit
# %%

# Load the calibration data
