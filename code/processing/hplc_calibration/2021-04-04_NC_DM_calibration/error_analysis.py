#%%
import numpy as np 
import pandas as pd
import altair as alt 
from altair_saver import save
import diaux.viz 
colors, palette = diaux.viz.altair_style()

DATA_PATH = '../../../../data/hplc_calibration/2021-04-05_NC_DM_calibration/processed'
data = pd.read_csv(f'{DATA_PATH}/2021-04-05_NC_DM_calibration_relative_areas.csv')
# %%
phos_data = data[(data['compound']=='phosphate') & (data['buffer_base']=='N-C-')]
avg_phos = phos_data['area'].mean()
phos_data['rel_area'] = phos_data['area'].values / avg_phos
phos_data['percent_err'] = (phos_data['area'] - avg_phos) / avg_phos
phos_data['sample_ID'] = np.arange(len(phos_data)) + 1
base = alt.Chart(phos_data, width=600, height=200).encode(
            x=alt.X('sample_ID:N', title='injection replicate')
)

area = base.mark_point(size=80).encode(
        y=alt.Y('area:Q', title='integrated signal [mV]',
            scale=alt.Scale(domain=[3E6, 6E6]))
)

rel_area = base.mark_point(size=80).encode(
        y=alt.Y('rel_area:Q', title='signal relative to mean')
)

err_area = base.mark_point(size=80).encode(
        y=alt.Y('percent_err:Q', title='percent error of mean',
                axis=alt.Axis(format='%')))



layer = area & rel_area & err_area
layer.properties(title='phosphate peak error analysis')
save(layer, './output/2021-04-05_NC_phosphate_calibration_error.pdf')
layer

# %%
