#%%
import numpy as np 
import pandas as pd
import scipy.stats
import diaux.viz
import altair as alt
from altair_saver import save
colors, palette = diaux.viz.altair_style()
alt.data_transformers.enable('json')

# Load the chromatograms for the calibration and diauxic shift data
cal_chroms = pd.read_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_calibration_chromatograms.csv')
diaux_chroms = pd.read_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_diauxic_shift_chromatograms.csv')

#%% Generate figures of the calibration curves and diauxic shift chromatograms
gluc_chrom = alt.Chart(cal_chroms[cal_chroms['carbon_source']=='glucose']
                       ).mark_line(size=1).encode(
                           x=alt.X('time_min:Q', title='retention time [min]'),
                           y=alt.Y('intensity_mV:Q', title='intensity [mV]'),
                           color=alt.Color('concentration_mM:O', title='glucose [mM]')
                       )

acetate_chrom = alt.Chart(cal_chroms[cal_chroms['carbon_source']=='acetate']
                      ).mark_line(size=1).encode(
                           x=alt.X('time_min:Q', title='retention time [min]'),
                           y=alt.Y('intensity_mV:Q', title='intensity [mV]'),
                           color=alt.Color('concentration_mM:O', title='acetate [mM]')
                       )

shift_chrom = alt.Chart(diaux_chroms).mark_line(size=1).encode(
                           x=alt.X('time_min:Q', title='retention time [min]'),
                           y=alt.Y('intensity_mV:Q', title='intensity [mV]'),
                           color=alt.Color('sample_id:N', title='sample ID')
                       )

layer = (gluc_chrom & acetate_chrom & shift_chrom).resolve_scale(color='independent')
save(layer, './output/2021-02-15_hplc_chromatograms.pdf')

#%%
# Perform background subtraction on a flat region of the chromatogram

# For the calibration data
dfs = []
for g, d in cal_chroms.groupby(['carbon_source', 'concentration_mM']):
    avg = d[(d['time_min'] >= 17) & (d['time_min']<=23)]['intensity_mV'].mean()
    d['corrected_intensity_mV']= d['intensity_mV'].values - avg
    dfs.append(d)
cal_chroms = pd.concat(dfs, sort=False)

# For the diauxic shift data
dfs = []
for g, d in diaux_chroms.groupby(['sample_id']):
    avg = d[(d['time_min'] >= 18) & (d['time_min']<=22)]['intensity_mV'].mean()
    d['corrected_intensity_mV']= d['intensity_mV'].values - avg
    dfs.append(d)
diaux_chroms = pd.concat(dfs, sort=False)

#%% Generate a plot for the peaks of interest
glucose_range = cal_chroms[(cal_chroms['carbon_source']=='glucose') & 
                           (cal_chroms['time_min'] >= 14) &
                           (cal_chroms['time_min'] <=17)]

acetate_range = cal_chroms[(cal_chroms['carbon_source']=='acetate') & 
                           (cal_chroms['time_min'] >= 24) &
                           (cal_chroms['time_min'] <=26)]
shift_range = diaux_chroms[(diaux_chroms['time_min'] >=14) & 
                            (diaux_chroms['time_min'] <=26)]
glucose = alt.Chart(glucose_range).mark_line().encode(
                    x=alt.X(field='time_min', type='quantitative', title='retention time [min]'),
                    y=alt.Y(field='corrected_intensity_mV', type='quantitative', title='subtracted intensity [mV]'),
                    color=alt.Color(field='concentration_mM', type='ordinal', title='glucose [mM]')
                    )
acetate = alt.Chart(acetate_range).mark_line().encode(
                    x=alt.X(field='time_min', type='quantitative', title='retention time [min]'),
                    y=alt.Y(field='corrected_intensity_mV', type='quantitative', title='subtracted intensity [mV]'),
                    color=alt.Color(field='concentration_mM', type='ordinal', title='acetate [mM]')
                    )
acetate = alt.Chart(acetate_range).mark_line().encode(
                    x=alt.X(field='time_min', type='quantitative', title='retention time [min]'),
                    y=alt.Y(field='corrected_intensity_mV', type='quantitative', title='subtracted intensity [mV]'),
                    color=alt.Color(field='concentration_mM', type='ordinal', title='acetate [mM]')
                    )
shift = alt.Chart(shift_range).mark_line().encode(
                    x=alt.X(field='time_min', type='quantitative', title='retention time [min]'),
                    y=alt.Y(field='corrected_intensity_mV', type='quantitative', title='subtracted intensity [mV]'),
                    color=alt.Color(field='sample_id', type='nominal', title='timepoint')
                    )
save((glucose | acetate | shift ).resolve_scale(color='independent'), 
    './output/2021-02-15_hplc_glucose_acetate_peaks.pdf') 

#%%
# Do my own integration for the peaks given this subtraction
bounds = {'glucose': [15.1, 16.1],
          'acetate': [24.3, 25.6]}
calibration = pd.DataFrame([])
for g, d in cal_chroms.groupby(['carbon_source', 'concentration_mM']):
    _d = d[(d['time_min']>=bounds[g[0]][0]) & (d['time_min']<=bounds[g[0]][1])]
    intensity = _d['corrected_intensity_mV'].sum()
    calibration = calibration.append({'integrated_intensity':intensity,
                                      'carbon_source':g[0],
                                      'concentration_mM':g[1],
                                      'arrival_time':bounds[g[0]][0],
                                      'departure_time':bounds[g[0]][1]},
                                      ignore_index=True)

# Integrate the peak areas for glucose and acetate for each shift sample
shift = pd.DataFrame([])
for idx, d in diaux_chroms.groupby(['sample_id']):
    for g in bounds.keys():
        _d = d[(d['time_min']>=bounds[g][0]) & (d['time_min']<=bounds[g][1])]
        intensity = _d['corrected_intensity_mV'].sum()
        shift = shift.append({'integrated_intensity':intensity,
                                      'carbon_source':g,
                                      'sample_id': idx,
                                      'arrival_time':bounds[g][0],
                                      'departure_time':bounds[g][1]},
                                      ignore_index=True)


#%%
# Do a simple linear regression to get a calibration curve
calibration_stats = pd.DataFrame([])
fit_dict = {}
fit_dfs = []
conc = np.linspace(0, 7, 200)
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

calibration.to_csv('./output/2021-02-15_acetate_glucose_calibration.csv', index=False) 
fit_df = pd.concat(fit_dfs, sort=False)

# Plot the calibration curve
cal_points = alt.Chart(calibration).mark_point(size=80).encode( 
                x=alt.X(field='concentration_mM', type='quantitative', title='concentration [mM]'),
                y=alt.Y(field='integrated_intensity', type='quantitative', title='peak intensity [mV]'),
                color=alt.Color(field='carbon_source', type='nominal', title='carbon source')
)
cal_fit = alt.Chart(fit_df).mark_line(clip=True).encode( 
                x=alt.X(field='concentration_mM', type='quantitative', title='concentration [mM]'),
                y=alt.Y(field='integrated_intensity', type='quantitative', title='peak intensity [mV]',
                        scale=alt.Scale(domain=[0, 3E5])),
                color=alt.Color(field='label', type='nominal', title='linear relationship',
                                scale=alt.Scale(domain=fit_df['label'].unique(),
                                                  range=palette[:2]))
)
layer = (cal_fit + cal_points).resolve_scale(color='independent')
layer
save(layer, './output/2021-02-15_calibration_curve.pdf')
#%%
# Given the calibration curve, compute the measured concentrations in the samples
dfs = []
for g, d in shift.groupby(['carbon_source', 'sample_id']):
    
    _integrated_intensity = d['integrated_intensity'].values[0]
    conc = (_integrated_intensity - fit_dict[g[0]]['intercept']) / fit_dict[g[0]]['slope']
    d['concentration_mM'] = conc
    dfs.append(d)

shift_conc = pd.concat(dfs, sort=False)

# Replace negative regions with zeros
shift_conc.loc[shift_conc['concentration_mM'] < 0, 'concentration_mM'] = 0
shift_conc['rescaled_concentration_mM'] = shift_conc['concentration_mM'].values * 5

# Save to csv
shift_conc.to_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_diauxic_shift_glucose_acetate_turnover.csv', index=False)

# %%
