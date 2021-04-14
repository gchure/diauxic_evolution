#%%
import numpy as np 
import cremerlab.hplc 
import diaux.hplc
import glob
import pandas as pd 
import matplotlib.pyplot as plt 
import cmdstanpy 
import arviz as az
import imp 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
imp.reload(cremerlab.hplc)
imp.reload(diaux.hplc)
import diaux.viz 
colors, palette = diaux.viz.matplotlib_style()
#%%
DATA_PATH = '../../../../data/hplc_calibration/2021-04-05_NC_DM_calibration'
raw_files = glob.glob(f'{DATA_PATH}/raw/*.txt')
cremerlab.hplc.convert(raw_files, peak_table=True)
#%%
dfs = []
for buff, base in zip(['NC', 'DM'], ['N-C-', 'DM']):
    peak_files = glob.glob(f'{DATA_PATH}/raw/converted/*{buff}*peak_table.csv')
    for p in peak_files:
        _, _, _, _, _, conc, _, _ = p.split('/')[-1].split('_')
        conc = float(conc.split('mM')[0])
        _df = pd.read_csv(p, comment="#")
        _df['carbon_conc_mM'] = conc
        _df['buffer_base'] = base
        dfs.append(_df)
peaks = pd.concat(dfs, sort=False)

#%%
assigned_peaks = diaux.hplc.assign_peaks(peaks)

# Compute relative areas
rel_df = []
compounds = ['phosphate', 'NaCl']
for g, d in assigned_peaks.groupby(['carbon_conc_mM', 'buffer_base']):
    for c in compounds:
        denom = d[d['compound']==c]['area'].values[0]
        d[f'rel_area_{c}'] = d['area'].values / denom
    rel_df.append(d)
rel_area_df = pd.concat(rel_df, sort=False)
rel_area_df.to_csv(f'{DATA_PATH}/processed/2021-04-05_NC_DM_calibration_relative_areas.csv',
    index=False)

#%%
