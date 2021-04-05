#%%
import pandas as pd
import cremerlab.hplc
import glob
import matplotlib.pyplot as plt
DATA_PATH = '../../../data/2021-04-03_HPLC_buffer_calibration'

# Get the file list and convert
raw_files = glob.glob(f'{DATA_PATH}/raw/*.txt')
cremerlab.hplc.convert(raw_files, output_dir=f'{DATA_PATH}/processed/')

#%%
# N-C- Buffer
NC_files = sorted(glob.glob(f'{DATA_PATH}/processed/*_NC_*.csv'))
NC_chroms, NC_peaks, plot = cremerlab.hplc.batch_process(NC_files, time_window=[10, 27],
                                                   show_viz=True, **dict(buffer=50))
for ext in ['pdf', 'png']:
    plt.savefig(f'output/2021-04-03_NC_buffer_chromatograms.{ext}')

#%% DM Buffer
DM_files = sorted(glob.glob(f'{DATA_PATH}/processed/*_DM_*.csv'))
DM_chroms, DM_peaks, plot = cremerlab.hplc.batch_process(DM_files, time_window=[13, 27],
                                                   show_viz=True, **dict(buffer=50))
for ext in ['pdf', 'png']:
    plt.savefig(f'output/2021-04-03_DM_buffer_chromatograms.{ext}')

# %%
# Assign peak identiies
peak_table = []
for peaks, buffer in zip([NC_peaks, DM_peaks], ['N-C-', 'DM']):
    peaks.loc[(peaks['retention_time'] >= 10.8) & 
              (peaks['retention_time'] <= 11.1), 'compound'] = 'NaCl'
    peaks.loc[(peaks['retention_time'] >= 15.0) & 
              (peaks['retention_time'] <= 15.3), 'compound'] = 'phosphate'
    peaks.loc[(peaks['retention_time'] >= 15.35) & 
              (peaks['retention_time'] <= 15.85), 'compound'] = 'glucose'
    peaks.loc[(peaks['retention_time'] >= 21.5) & 
              (peaks['retention_time'] <= 22.1), 'compound'] = 'glycerol'
    peaks.loc[(peaks['retention_time'] >= 24.5) & 
              (peaks['retention_time'] <= 25.1), 'compound'] = 'acetate'

    # Parse the sample names and extract concentration
    peaks['carbon_conc_mM']  = [float(s.split('_')[-1].split('mM')[0]) for s in peaks['sample'].values]
    peaks['buffer_base'] = buffer 
    peak_table.append(peaks)
peaks = pd.concat(peak_table, sort=False)

# Drop unnecessary columns and save
peaks.drop(columns=['sample'], inplace=True)
# peaks.dropna(inplace=True)
peaks.to_csv(f'{DATA_PATH}/processed/2021-04-03_calibration_peaks.csv', index=False)

#%%
# Compute the relative  areas
rel_peaks = []
for g, d in peaks.groupby(['buffer_base', 'carbon_conc_mM']):
    phos = d[d['compound']=='phosphate']
    unk = d[d['compound']=='NaCl']
    if len(unk) != 0:
        iterator = zip([phos['area'].values[0], unk['area'].values[0]], 
                        ['phosphate', 'NaCl'])
    else:
        iterator = zip([phos['area'].values[0]], ['phosphate'])
    for p, n in iterator:
        d = d.copy(deep=True)
        d['rel_area'] = d['area'].values / p
        d['relative_peak'] = n
        rel_peaks.append(d)

rel_peaks = pd.concat(rel_peaks, sort=False)
rel_peaks = rel_peaks[rel_peaks['compound'].isin(
                    ['glucose', 'glycerol', 'acetate'])]
# %%
rel_peaks.to_csv(f'{DATA_PATH}/processed/2021-04-03_relative_calibration_peaks.csv', 
                 index=False)

# %%
