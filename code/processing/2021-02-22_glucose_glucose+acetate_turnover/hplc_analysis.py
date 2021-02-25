#%%
import numpy as np 
import pandas as pd
import scipy.stats
import diaux.viz
import altair as alt
from altair_saver import save
colors, palette = diaux.viz.altair_style()
alt.data_transformers.enable('json')
DATA_PATH = '../../../data'
DATE = '2021-02-22'
STRAIN = 'NCM3722'
CONDITION = 'glucose_glucose+acetate_turnover'
# Load the chromatograms for the calibration and diauxic shift data
cal_chroms = pd.read_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}/processed/{DATE}_{STRAIN}_calibration_chromatograms.csv')
turnover_chroms = pd.read_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}/processed/{DATE}_{STRAIN}_{CONDITION}_chromatograms.csv')

#%%
glucose_range = cal_chroms[(cal_chroms['time_min'] >= 14.5) &
                           (cal_chroms['time_min'] <=16.5)]

acetate_range = cal_chroms[(cal_chroms['time_min'] >= 24) &
                           (cal_chroms['time_min'] <=26)]

# Plot the raw chromatogram chroms
cal_tot = alt.Chart(cal_chroms, width=600).mark_line(size=1).encode(
                    x=alt.X('time_min:Q', title='retention time [min]'),
                    y=alt.Y('intensity_mV:Q', title='intensity [mV]'),
                    color=alt.Color('glucose_mM:O', legend=None)
                    ).properties(title='calibration chromatograms')
gluc_tot = alt.Chart(glucose_range, width=200, height=200).mark_line(size=1).encode(
                    x=alt.X('time_min:Q', title='retention time [min]'),
                    y=alt.Y('intensity_mV:Q', title='intensity [mV]'),
                    color=alt.Color('glucose_mM:O', title='glucose [mM]')
                    ).properties(title='glucose titration')
ace_tot = alt.Chart(acetate_range, width=200, height=200).mark_line(size=1).encode(
                    x=alt.X('time_min:Q', title='retention time [min]'),
                    y=alt.Y('intensity_mV:Q', title='intensity [mV]'),
                    color=alt.Color('acetate_mM:O', title='acetate [mM]')
                    ).properties(title='acetate titration')
layer = (cal_tot & (gluc_tot | ace_tot).resolve_scale(color='independent')).resolve_scale(color='independent')
save(layer, f'./output/{DATE}_calibration_chromatograms.pdf')
layer
#%%
# Perform background subtraction on a flat region of the chromatogram

# For the calibration data
dfs = []
for g, d in cal_chroms.groupby(['glucose_mM']):
    avg = d[(d['time_min'] >= 17) & (d['time_min']<=23)]['intensity_mV'].mean()
    d['corrected_intensity_mV']= d['intensity_mV'].values - avg
    dfs.append(d)
cal_chroms = pd.concat(dfs, sort=False)

# For the diauxic shift data
dfs = []
for g, d in turnover_chroms.groupby(['replicate', 'sample_id', 'glucose_mM', 'acetate_mM']):
    avg = d[(d['time_min'] >= 18) & (d['time_min']<=22)]['intensity_mV'].mean()
    d['corrected_intensity_mV']= d['intensity_mV'].values - avg
    dfs.append(d)
turnover_chroms = pd.concat(dfs, sort=False)

# Plot the subtracted 
glucose_range = cal_chroms[(cal_chroms['time_min'] >= 14.5) &
                           (cal_chroms['time_min'] <=16.5)]

acetate_range = cal_chroms[(cal_chroms['time_min'] >= 24) &
                           (cal_chroms['time_min'] <=26)]

# Plot the raw chromatogram chroms
cal_tot = alt.Chart(cal_chroms, width=600).mark_line(size=1).encode(
                    x=alt.X('time_min:Q', title='retention time [min]'),
                    y=alt.Y('corrected_intensity_mV:Q', title='corrected intensity [mV]'),
                    color=alt.Color('glucose_mM:O', legend=None)
                    ).properties(title='calbiration chromatograms')
gluc_tot = alt.Chart(glucose_range, width=200, height=200).mark_line(size=1).encode(
                    x=alt.X('time_min:Q', title='retention time [min]'),
                    y=alt.Y('corrected_intensity_mV:Q', title='corrected intensity [mV]'),
                    color=alt.Color('glucose_mM:O', title='glucose [mM]')
                    ).properties(title='glucose titration')
ace_tot = alt.Chart(acetate_range, width=200, height=200).mark_line(size=1).encode(
                    x=alt.X('time_min:Q', title='retention time [min]'),
                    y=alt.Y('corrected_intensity_mV:Q', title='corrected intensity [mV]'),
                    color=alt.Color('acetate_mM:O', title='acetate [mM]')
                    ).properties(title='acetate titration')
layer = (cal_tot & (gluc_tot | ace_tot).resolve_scale(color='independent')).resolve_scale(color='independent')
save(layer, f'./output/{DATE}_calibration_chromatograms_subtracted.pdf')
layer
#%%
# Do my own integration for the peaks given this subtraction
bounds = {'glucose': [15.1, 16.1],
          'acetate': [24.3, 25.6]}
calibration = pd.DataFrame([])
for g, d in cal_chroms.groupby(['glucose_mM', 'acetate_mM']):
    for k, v in bounds.items(): 
        _d = d[(d['time_min']>=v[0]) & (d['time_min']<=v[1])]
        intensity = _d['corrected_intensity_mV'].sum()
        if k == 'glucose':
            conc = g[0]
        else:
            conc = g[1]
        calibration = calibration.append({'integrated_intensity':intensity,
                                      'carbon_source':k,
                                      'concentration_mM':conc,
                                      'arrival_time':v[0],
                                      'departure_time':v[1]},
                                      ignore_index=True)

# Integrate the peak areas for glucose and acetate for each turnover sample
turnover = pd.DataFrame([])
for g, d in turnover_chroms.groupby(['glucose_mM', 'acetate_mM', 'replicate', 'sample_id']):
    for k, v  in bounds.items():
        _d = d[(d['time_min']>=v[0]) & (d['time_min']<=v[1])]
        intensity = _d['corrected_intensity_mV'].sum()
        turnover = turnover.append({'integrated_intensity':intensity,
                                      'carbon_source':k,
                                      'replicate': g[-2],
                                      'sample_id': g[-1],
                                      'arrival_time':v[0],
                                      'departure_time':v[1],
                                      'glucose_mM': g[0],
                                      'acetate_mM': g[1]},
                                      ignore_index=True)


#%%
# Do a simple linear regression to get a calibration curve
calibration_stats = pd.DataFrame([])
fit_dict = {}
fit_dfs = []
conc = np.linspace(0, 30, 200)
for g, d in calibration.groupby(['carbon_source']):
    popt = scipy.stats.linregress(d['concentration_mM'], d['integrated_intensity'])
    fit = popt[0] * conc + popt[1]
    _df = pd.DataFrame([])
    _df['concentration_mM'] = conc
    _df['integrated_intensity'] = fit
    _df['carbon_source'] = g
    fit_dict[g] = {'slope':popt[0], 'intercept':popt[1]}  
    _df['label'] = f'{popt[0]:0.0f}[c] + {popt[1]:0.0f}'
    fit_dfs.append(_df) 

    # Add the inference data to the calibration
    calibration.loc[calibration['carbon_source']==g, 'slope_mV_mM'] = popt[0]
    calibration.loc[calibration['carbon_source']==g, 'intercept_mV_mM'] = popt[1]

calibration.to_csv(f'./output/{DATE}_acetate_glucose_calibration.csv', index=False) 
fit_df = pd.concat(fit_dfs, sort=False)

# Plot the calibration curve
cal_points = alt.Chart(calibration).mark_point(size=80).encode( 
                x=alt.X(field='concentration_mM', type='quantitative', title='concentration [mM]'),
                y=alt.Y(field='integrated_intensity', type='quantitative', title='peak intensity [mV]'),
                color=alt.Color(field='carbon_source', type='nominal', title='carbon source')
)
cal_fit = alt.Chart(fit_df).mark_line(clip=True).encode( 
                x=alt.X(field='concentration_mM', type='quantitative', title='concentration [mM]'),
                y=alt.Y(field='integrated_intensity', type='quantitative', title='peak intensity [mV]'),
                        
                color=alt.Color(field='label', type='nominal', title='linear relationship',
                                scale=alt.Scale(domain=fit_df['label'].unique(),
                                                  range=palette[:2]))
)
layer = (cal_fit + cal_points).resolve_scale(color='independent')
layer
save(layer, f'./output/{DATE}_calibration_curve.pdf')
#%%
# # Given the calibration curve, compute the measured concentrations in the samples
# dfs = []
# for g, d in shift.groupby(['carbon_source', 'sample_id']):
    
#     _integrated_intensity = d['integrated_intensity'].values[0]
#     conc = (_integrated_intensity - fit_dict[g[0]]['intercept']) / fit_dict[g[0]]['slope']
#     d['concentration_mM'] = conc
#     dfs.append(d)

# shift_conc = pd.concat(dfs, sort=False)

# # Replace negative regions with zeros
# shift_conc.loc[shift_conc['concentration_mM'] < 0, 'concentration_mM'] = 0
# shift_conc['rescaled_concentration_mM'] = shift_conc['concentration_mM'].values * 5

# # Save to csv
# shift_conc.to_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_diauxic_shift_glucose_acetate_turnover.csv', index=False)

# # %%
