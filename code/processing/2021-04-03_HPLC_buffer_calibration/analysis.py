#%%
import sys 
sys.path.insert(0, '../../../')
import numpy as np 
import pandas as pd 
import diaux.viz 
import altair as alt
colors, palette = diaux.viz.altair_style()
# %%
# Load the two calibration datasets
DATA_PATH = '../../../data/2021-04-03_HPLC_buffer_calibration/processed'
peaks = pd.read_csv(f'{DATA_PATH}/2021-04-03_calibration_peaks.csv')
rel_peaks = pd.read_csv(f'{DATA_PATH}/2021-04-03_relative_calibration_peaks.csv')

# %%
peak_base = alt.Chart(peaks).encode(
                x=alt.X('carbon_conc_mM:Q', title='carbon source concentration [mM]'),
                y=alt.Y('area:Q', title='integrated signal [mV]'),
                shape=alt.Shape('buffer_base:N', title='buffer base'),
                color=alt.Color('compound:N', title='compound')
)

(peak_base.mark_point(size=80, opacity=0.5))

# %%
# Look first at N-C-
nc = rel_peaks[rel_peaks['buffer_base']=='N-C-']
nacl = nc[nc['relative_peak']=='NaCl']
nacl

# %%
alt.Chart(nacl).mark_point().encode(
        x=alt.X('carbon_conc_mM:Q', title='carbon source concentration [mM]'),
        y=alt.Y('rel_area:Q', title='signal relative to NaCl'),
        color=alt.Color('compound:N')
)

# %%
