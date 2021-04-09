#%%
import numpy as np 
import cremerlab.hplc 
import glob
import pandas as pd 
import matplotlib.pyplot as plt 
import imp 
imp.reload(cremerlab.hplc)

DATA_PATH = '../../../../data/hplc_calibration/2021-04-05_NC_DM_calibration'
raw_files = glob.glob(f'{DATA_PATH}/raw/*.txt')
cremerlab.hplc.convert(raw_files, peak_table=True)
#%%
peak_files = glob.glob(f'{DATA_PATH}/raw/converted/*peak_table.csv')
peak_df = []
for p in peak_files:
    _, buffer_base, _, _, _, conc, _, _ = p.split('/')[-1].split('_')
    conc = float(conc.split('mM')[0])
    _df = pd.read_csv(p, comment="#")
    _df['carbon_conc_mM'] = conc
    _df['buffer_base'] = buffer_base
    peak_df.append(_df)
peaks = pd.concat(peak_df, sort=False)

#%%
# Assign the peaks
peaks.loc[(peaks['retention_time'] >=10) & (peaks['retention_time'] <= 11.5),
            'compound'] = 'NaCl'
peaks.loc[(peaks['retention_time'] >=13.5) & (peaks['retention_time'] <= 14.1),
            'compound'] = 'phosphate'

peaks.loc[(peaks['retention_time'] >=24) & (peaks['retention_time'] <= 26),
            'compound'] = 'acetate'
peaks.loc[(peaks['retention_time'] >=21.9) & (peaks['retention_time'] <= 22.3),
            'compound'] = 'glycerol'
peaks.loc[(peaks['retention_time'] >=15) & (peaks['retention_time'] <= 16),
            'compound'] = 'glucose'
peaks.to_csv(f'{DATA_PATH}/processed/2021-04-05_NC_DM_carbon_source_calibration.csv', 
              index=False)
#%%
files = sorted(glob.glob(f'{DATA_PATH}/raw/converted/*NC*chromatogram.csv'))
chroms, peaks, plot = cremerlab.hplc.batch_process(files, time_window=[10, 27],
                                                   show_viz=True, **dict(buffer=50))
plt.savefig('./output/2021-04-4_NC_calibration_chromatograms.pdf')
# %%
for i, p in enumerate([peaks, peak_df]):
    p.loc[(p['retention_time'] >= 10.5) & (p['retention_time'] <=11.1),
               'compound'] = 'NaCl'
    p.loc[(p['retention_time'] >= 14.0) & (p['retention_time'] <=15.0),
               'compound'] = 'Phosphate'
    p.loc[(p['retention_time'] >= 14.0) & (p['retention_time'] <=15.0),
               'compound'] = 'Phosphate'   
    p.loc[(p['retention_time'] >=23) & (p['retention_time'] <= 25),
                'compound'] = 'acetate'
    p.loc[(p['retention_time'] >=20) & (p['retention_time'] <= 22.5),
                'compound'] = 'glycerol'
    p.loc[(p['retention_time'] >= 15.3) & (p['retention_time'] <= 16),
                'compound'] = 'glucose'

    if i == 0: 
        p['carbon_conc'] = [float(s.split('_')[-2].split('mM')[0]) for s in peaks['sample'].values]
    p.dropna(inplace=True)

# %%
import sys
sys.path.insert(0, '../../../../')
import altair as alt 
from altair_saver import save 
import diaux.viz
colors = diaux.viz.altair_style()

#%%
chart = alt.Chart(peak_df).mark_point(size=80, opacity=0.75).encode(
            x=alt.X('carbon_conc_mM:Q', title='carbon source concentration [mM]'),
            y=alt.Y('area:Q', title='integrated signal [mV]'),
            color=alt.Color('compound:N', title='compound ID')
)

chart

#%%
rel_table = []
for g, d in peak_df.groupby(['carbon_conc_mM']):
    nacl = d[d['compound']=='Phosphate']['area'].values[0]
    d['rel_area'] = d['area'].values / nacl
    rel_table.append(d)
rel_table = pd.concat(rel_table, sort=False)
chart = alt.Chart(rel_table).mark_point(size=80, opacity=0.75).encode(
            x=alt.X('carbon_conc_mM:Q', title='carbon source concentration [mM]'),
            y=alt.Y('rel_area:Q', title='relative signal'),
            color=alt.Color('compound:N', title='compound ID')
)

chart
# %%
chart = alt.Chart(peaks).mark_point(size=80, opacity=0.75).encode(
            x=alt.X('carbon_conc:Q', title='carbon source concentration [mM]'),
            y=alt.Y('area:Q', title='integrated signal [mV]'),
            color=alt.Color('compound:N', title='compound ID')
)

save(chart, './output/2021-04-04_NC_calibration.pdf')
# %%
rel_df = []
for g, d in peaks.groupby(['carbon_conc']):
    nacl = d[d['compound']=='NaCl']['area'].values[0]
    d['rel_area'] = d['area'].values / nacl
    rel_df.append(d)
rel_peaks = pd.concat(rel_df, sort=False)
# %%
chart = alt.Chart(rel_peaks).mark_point(size=80, opacity=0.75).encode(
            x=alt.X('carbon_conc:Q', title='carbon source concentration [mM]'),
            y=alt.Y('rel_area:Q', title='relative signal'),
            color=alt.Color('compound:N', title='compound ID')
)

save(chart, './output/2021-04-04_NC_rel_calibration.pdf')
# %%
