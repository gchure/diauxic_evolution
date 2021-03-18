#%%
import numpy as np 
import pandas as pd 
import altair as alt 
import cmdstanpy
import arviz as az
from altair_saver import save
import diaux.viz
colors, palette = diaux.viz.altair_style()
# %%
# Load the peak table
peaks =  pd.read_csv('../../../data/2021-03-15_HPLC_calibration/processed/2021-03-15_calibration_peaks.csv')
peaks.dropna(inplace=True)

peaks['retention_time'] = peaks['retention_time'].values.astype(float)
peaks['initial_time'] = peaks['initial_time'].values.astype(float)
peaks['final_time'] = peaks['final_time'].values.astype(float)

# Restrict the retention times between 10 m and 30
peaks = peaks[(peaks['retention_time'] >= 10) & (peaks['retention_time'] <= 30)]

# Ignore things below a milimolar
peaks = peaks[(peaks['solute_concentration_mM'] >= 0.5) & (peaks['solute_concentration_mM'] < 30)]

# Categorize based on retention time.
peaks.loc[(peaks['retention_time'].round(decimals=0) >=15) & 
          (peaks['retention_time'].round(decimals=0) <= 16), 'solute'] = 'glucose'
peaks.loc[(peaks['retention_time'].round(decimals=0) >= 22) & 
          (peaks['retention_time'].round(decimals=0) <= 23), 'solute'] = 'glycerol'
peaks.loc[(peaks['retention_time'].round(decimals=0) >= 24) & 
          (peaks['retention_time'].round(decimals=0) <= 25), 'solute'] = 'acetate'

# Remove peaks that couldn't be assigned
peaks.dropna(inplace=True)

#%%
samples_dfs, stats_dfs = [], []

# Load the inferrential model
model = cmdstanpy.CmdStanModel(stan_file='../../stan/general_linear_model.stan')
quantiles = [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
for g, d in peaks.groupby('solute'):
    data_dict = {'N': len(d),
                 'x':d['solute_concentration_mM'].values.astype(float),
                 'y':d['peak_area'].values.astype(float)}

    # Sample and convert to a dataframe
    samples = model.sample(data=data_dict)
    samples = az.from_cmdstanpy(samples)
    samples = samples.posterior.to_dataframe().reset_index()
    samples = samples[['slope', 'intercept', 'sigma']]
    samples['solute'] = g
    samples_dfs.append(samples)

    # Compute the summary
    samp_quants = samples.quantile(quantiles).reset_index()
    samp_quants.rename(columns={'index':'quantile'}, inplace=True)
    samp_quants['solute'] = g
    stats_dfs.append(samp_quants)

#%%
samples = pd.concat(samples_dfs, sort=False)
stats = pd.concat(stats_dfs, sort=False )
samples.to_csv('../../../data/2021-03-15_HPLC_calibration/processed/2021-03-15_HPLC_calibration_mcmc_samples.csv')
stats.to_csv('../../../data/2021-03-15_HPLC_calibration/processed/2021-03-15_HPLC_calibration_mcmc_statistics.csv')
#%%
# Compute fits. 
conc_range = np.linspace(0, 30, 100)
fits = []
for g, d in samples.groupby(['solute']): 
    cred_region = np.zeros((2, len(conc_range)))
    medians = []
    _samples = samples[samples['solute']==g]
    for i, c in enumerate(conc_range):
        _fit = _samples['slope'].values * c + _samples['intercept'].values
        median =np.median(_fit) 
        cred_region[:, i] = np.percentile(_fit, [2.5, 97.5])
        medians.append(median)
    _df = pd.DataFrame([])
    _df['solute_concentration_mM'] = conc_range
    _df['lower'] = cred_region[0, :]
    _df['upper'] = cred_region[1, :]
    _df['median'] = medians
    _df['solute'] = g
    fits.append(_df)
fits = pd.concat(fits, sort=False)

#%% Set up the plot
glucose_data = peaks[peaks['solute']=='glucose']
glucose_fit = fits[fits['solute']=='glucose']
glycerol_data = peaks[peaks['solute']=='glycerol']
glycerol_fit = fits[fits['solute']=='glycerol']
acetate_data = peaks[peaks['solute']=='acetate']
acetate_fit = fits[fits['solute']=='acetate']

#%%
ac_stats = stats[(stats['solute']=='acetate')]
glc_stats = stats[(stats['solute']=='glucose')]
gly_stats = stats[(stats['solute']=='glycerol')]
acetate_title = f"acetate; α ≈ {int(ac_stats[ac_stats['quantile']==0.50]['slope'].values[0].round(-2))} mV / mM"
glucose_title = f"glucose; α ≈ {int(glc_stats[glc_stats['quantile']==0.50]['slope'].values[0].round(-2))} mV / mM"
glycerol_title = f"glycerol; α ≈ {int(gly_stats[gly_stats['quantile']==0.50]['slope'].values[0].round(-2))} mV / mM"

#%%
glucose_points = alt.Chart(glucose_data).mark_point(size=80, color=colors['blue']).encode(
                    x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
                    y=alt.Y('peak_area:Q', title='peak area [mV]'))
glucose_lines = alt.Chart(glucose_fit
                          ).mark_area(
                               color=colors['blue'], opacity=0.25
                          ).encode(
                    x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
                    y=alt.Y('lower:Q', title='peak area [mV]'),
                    y2='upper:Q')
glucose_median = alt.Chart(glucose_fit
                          ).mark_line(
                               color=colors['blue']
                          ).encode(
                    x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
                    y=alt.Y('median:Q', title='peak area [mV]'))


glycerol_points = alt.Chart(glycerol_data).mark_point(size=80, color=colors['green']).encode(
                    x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
                    y=alt.Y('peak_area:Q', title='peak area [mV]'))
glycerol_lines = alt.Chart(glycerol_fit
                          ).mark_area(
                               color=colors['green'], opacity=0.25
                          ).encode(
                    x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
                    y=alt.Y('lower:Q', title='peak area [mV]'),
                    y2='upper:Q')
glycerol_median = alt.Chart(glycerol_fit
                          ).mark_line(
                               color=colors['green']
                          ).encode(
                    x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
                    y=alt.Y('median:Q', title='peak area [mV]'))

acetate_points = alt.Chart(acetate_data).mark_point(size=80, color=colors['black']).encode(
                    x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
                    y=alt.Y('peak_area:Q', title='peak area [mV]'))
acetate_lines = alt.Chart(acetate_fit
                          ).mark_area(
                               color=colors['black'], opacity=0.25
                          ).encode(
                    x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
                    y=alt.Y('lower:Q', title='peak area [mV]'),
                    y2='upper:Q')
acetate_median = alt.Chart(acetate_fit
                          ).mark_line(
                               color=colors['black']
                          ).encode(
                    x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
                    y=alt.Y('median:Q', title='peak area [mV]'))

glucose_plots = (glucose_lines + glucose_points + glucose_median).properties(title=glucose_title)
glycerol_plots = (glycerol_lines + glycerol_points + glycerol_median).properties(title=glycerol_title)
acetate_plots = (acetate_lines + acetate_points + acetate_median).properties(title=acetate_title)

out = (glucose_plots | glycerol_plots | acetate_plots)
save(out, './output/2021-03-15_calibration_plots.pdf')
# %%
alt.Chart(peaks).mark_point().encode(
        x=alt.X('solute_concentration_mM:Q', title='concentration [mM]'),
        y=alt.Y('peak_area:Q', title='peak area [mV]')  ,
        color=alt.Color('solute:N', title='solute')
)

# %%
