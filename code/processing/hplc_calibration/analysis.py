#%%
import pandas as pd 
import cmdstanpy 
import arviz as az

DATA_PATH = '../../../../data/hplc_calibration/2021-04-05_NC_DM_calibration'
rel_area_df = pd.read_csv(f'{DATA_PATH}/processed/2021-04-05_NC_DM_calibration_relative_areas.csv')
model = cmdstanpy.CmdStanModel(stan_file='../../../stan/general_linear_model.stan')
samples_df = []
for g, d in rel_area_df[rel_area_df['compound'].isin(
                        ['glucose', 'glycerol', 'acetate'])].groupby(
                        ['buffer_base', 'compound']):
    data_dict = {'N':len(d),
                 'x':d['carbon_conc_mM'].values.astype(float),
                 }
    for key in ['rel_area_phosphate', 'rel_area_NaCl', 'area']:
        data_dict['y'] = d[key].values.astype(float)
        samples = model.sample(data=data_dict)
        samples = az.from_cmdstanpy(samples)
        samples = samples.posterior.to_dataframe().reset_index()
        samples = samples[['slope', 'intercept', 'sigma']]
        samples['buffer_base'] = g[0]
        samples['compound'] = g[1]
        samples['area_param'] = key
        samples_df.append(samples)
samples_df = pd.concat(samples_df, sort=False)
samples_df.to_csv(f'{DATA_PATH}/processed/2021-04-05_NC_DM_calibration_MCMC.csv', index=False)



# %%
