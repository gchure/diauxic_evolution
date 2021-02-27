#%%
import numpy as np 
import pandas as pd 
import diaux.viz 
import arviz as az
import cmdstanpy
import altair as alt
from altair_saver import save
colors, _ = diaux.viz.altair_style()
alt.data_transformers.disable_max_rows()

# Define  constants of the experiment 
DATE = '2021-02-25'
STRAIN = 'NCM3722'
CONDITION = 'glucose'
OTHER = 'RNA_isolation'
DATA_PATH = '../../../data'

# Load the raw data 
data = pd.read_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}_{OTHER}/raw/{DATE}_{STRAIN}_{CONDITION}_growth.csv')

# Do minor datacleaning
data = data[(data['od_600nm'] >= 0.04) & (data['od_600nm'] <= 0.4)]

# Compute the elapsed time
data['clock_time'] = pd.to_datetime(data['clock_time'])
data.sort_values(by='clock_time', inplace=True) 
data['elapsed_time_hr'] = data['clock_time'] - data['clock_time'].values[0]
data['rel_od_600nm'] = data['od_600nm'].values / data['od_600nm'].values[0]
data['elapsed_time_hr'] = (data['elapsed_time_hr'].astype('timedelta64[m]'))/60
data.drop(columns=['clock_time'], inplace=True)

data.head()

# %%
# Load the stan model 
sm = cmdstanpy.CmdStanModel(stan_file='../../stan/growth_rate_glm.stan')

# Define the data dictionary
data_dict = {'N': len(data),
             'absorbance': data['od_600nm'].values,
             'time': data['elapsed_time_hr'].values}

#%% Sample the model
samples = sm.sample(data=data_dict)

# Translate samples to an arviz object for inspection
samples = az.from_cmdstanpy(posterior=samples)
samples_df = samples.posterior.to_dataframe().reset_index()

#%%
# Save the sampling statistics and the processed growth data
samples_df.to_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}_{OTHER}/processed/{DATE}_{STRAIN}_{CONDITION}_parameter_samples.csv', index=False)
data.to_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}_{OTHER}/processed/{DATE}_{STRAIN}_{CONDITION}_growth.csv', index=False)
#%%
# Compute all of the best fits for display, ignoring the homoscedastic error for now. 
time_range = np.linspace(0, 2.5, 2)
fit_df = pd.DataFrame()
for t in time_range:
    fit = samples_df['absorbance_0'] * np.exp(samples_df['lambda'] * t)
    lower, upper = np.percentile(fit, [2.5, 97.5])
    fit_df = fit_df.append({'elapsed_time_hr': t,
                            'od_600nm_lower': lower,
                            'od_600nm_upper': upper},
                            ignore_index=True)
# %%
# Plot the growth data 
data_points = alt.Chart(data, width=600, height=200).mark_point(size=100).encode(
                x=alt.X('elapsed_time_hr:Q', title='elapsed time [inv. hr]',
                        scale=alt.Scale(domain=[0, 2.6])),
                y=alt.Y('od_600nm:Q', title='optical density [a.u.]',
                         scale=alt.Scale(type='log', domain=[0.04, 0.4]))
)

fit = alt.Chart(fit_df).mark_area(opacity=0.5).encode(
                    x=alt.X('elapsed_time_hr:Q', title='elapsed time [inv. hr]'),
                    y=alt.Y('od_600nm_lower:Q', title='optical density [a.u.]',
                            scale=alt.Scale(type='log', domain=[0.04, 0.4])),
                    y2='od_600nm_upper:Q',
                )
growth_curve = (fit + data_points)

# Plot the parameter estimates
lambda_posterior = alt.Chart(samples_df, height=250, width=275).mark_bar().encode(
                    x=alt.X('lambda:Q', title='growth rate [inv. hr]', 
                            bin=alt.Bin(maxbins=100)),
                    y=alt.Y('count()', title='number of samples')
)

abs0_posterior = alt.Chart(samples_df, height=250, width=275).mark_bar().encode(
                    x=alt.X('absorbance_0:Q', title='initial absorbance [a.u.]', 
                            bin=alt.Bin(maxbins=100)),
                    y=alt.Y('count()', title='number of samples')
)
posteriors = (lambda_posterior | abs0_posterior)
layout = growth_curve & posteriors
save(layout,f'./output/{DATE}_{STRAIN}_growth_statistics.png')
save(layout,f'./output/{DATE}_{STRAIN}_growth_statistics.pdf')

# %%
