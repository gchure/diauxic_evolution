#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import imp
import diaux.viz 
imp.reload(diaux.viz)
colors, palette = diaux.viz.matplotlib_style()
DATA_PATH = '../../../../data/hplc_calibration/2021-04-05_NC_DM_calibration/processed'

# Load the calibration data nd the MCMC samples
data = pd.read_csv(f'{DATA_PATH}/2021-04-05_NC_DM_calibration_relative_areas.csv')
mcmc = pd.read_csv(f'{DATA_PATH}/2021-04-05_NC_DM_calibration_MCMC.csv')

# %%
# Set up the outputs
output_files = ['2021-04-05_DM_calibration.pdf', '2021-04-05_NC_calibration.pdf']
_colors = {'glucose':colors['primary_blue'],
           'glycerol': colors['primary_green'],
           'acetate': colors['primary_red']}
conc_range = np.linspace(0, 35, 100)
for o in output_files:
    with PdfPages(f'./output/{o}') as  pdf:
        for g, d in data[data['compound'].isin(['glucose', 'glycerol', 'acetate'])].groupby('compound'):
            fig, ax = plt.subplots(3, 2, figsize=(6,6), sharex=True)
            for i in range(2):
                ax[0, i].set_ylabel('integrated signal [mV]')
                ax[1, i].set_ylabel('signal relative to phosphate')
                ax[2, i].set_ylabel('signal relative to NaCl')

            counter = 0 

            for _g, _d in d.groupby('buffer_base'):
                # Compute the calibration envelope
                ax[0, counter].set_title(f'{_g} buffer')
                for j, key in enumerate(['area', 'rel_area_phosphate', 'rel_area_NaCl']):
                    _mcmc = mcmc[(mcmc['buffer_base']==_g) & (mcmc['compound']==g) &
                                (mcmc['area_param']==key)]
                    cred_region = np.zeros((2, len(conc_range)))

                    for k, c in enumerate(conc_range):
                        fit = _mcmc['slope'].values * c + _mcmc['intercept'].values
                        cred_region[:, k] = np.percentile(fit, [2.5, 97.5])


                    # Plot the data
                    ax[j, counter].plot(_d['carbon_conc_mM'], _d[key], 'o', color=_colors[g]) 
                    ax[j, counter].fill_between(conc_range, cred_region[0, :], cred_region[1, :], color=_colors[g], alpha=0.4)


                counter += 1
            ax[-1,0].set_xlabel('concentration [mM]')
            ax[-1,1].set_xlabel('concentration [mM]')
            fig.suptitle(f'{g} calibration curve')
            pdf.savefig(fig)
# %%
