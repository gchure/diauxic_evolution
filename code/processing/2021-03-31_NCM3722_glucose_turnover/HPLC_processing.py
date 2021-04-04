#%%
import matplotlib.pyplot as plt
import glob
import cremerlab.hplc
#%%
# Define the data path  and convert the files
DATA_PATH = '../../../data'
raw_files = glob.glob(f'{DATA_PATH}/2021-03-31_NCM3722_glucose_turnover/raw/*.txt')
raw_cal_files = glob.glob(f'{DATA_PATH}/2021-04-02_NC_buffer_glucose_glycerol_acetate_calibration/raw/*.txt')
cremerlab.hplc.convert(raw_files)
cremerlab.hplc.convert(raw_cal_files)
# %%
# Separately load the calibration and turnover experiments and batch process
cal_files = sorted(glob.glob(f'{DATA_PATH}/2021-04-02_NC_buffer_glucose_glycerol_acetate_calibration/raw/converted/*calibration*.csv'))
cal_chroms, cal_peaks, plot = cremerlab.hplc.batch_process(cal_files, 
                                                           time_window=[10, 27],
                                                           show_viz=True,
                                                           **dict(prominence=5E-4))
fig, ax = plot
for ext in ['pdf', 'png']:
    plt.savefig(f'./output/2021-04-02_calibration_chromatograms.{ext}', bbox_inches='tight')

#%%
# Tidy the calibration peak table to map glucose and acetate concentrations
cal_peaks.loc[(cal_peaks['retention_time'] >= 15.4) &
              (cal_peaks['retention_time'] <= 15.6), 'compound'] = 'glucose'
cal_peaks.loc[(cal_peaks['retention_time'] >= 24.5) &
              (cal_peaks['retention_time'] <= 30.1), 'compound'] = 'acetate'
cal_peaks.dropna(inplace=True)

# Remap concentrations as a column
concs = [float(s.split('_')[-1].split('mM')[0]) for s in cal_peaks['sample'].values]
cal_peaks['concentration_mM'] = concs
cal_peaks.drop(columns=['sample'], inplace=True)
cal_peaks.to_csv(f'{DATA_PATH}/processed/2021-03-31_glucose_acetate_calibration.csv', 
                 index=False)

# %%
# Process the turnover data
samp_files = sorted(glob.glob(f'{DATA_PATH}/raw/converted/*NCglucose*.csv'))
samp_chroms, samp_peaks, samp_plot = cremerlab.hplc.batch_process(samp_files,
                                            time_window=[13, 27], show_viz=True)
fig, ax = plot
for ext in ['pdf', 'png']:
    plt.savefig(f'./output/2021-03-31_turnover_chromatograms.{ext}', bbox_inches='tight')

# %%
# Add replicate and timepoint info to the peak table
repls = [int(s.split('_')[4]) for s in samp_peaks['sample'].values]
samp_idx = [int(s.split('_')[-1]) for s in samp_peaks['sample'].values]
samp_peaks['replicate'] = repls
samp_peaks['sample_idx'] = samp_idx
#%%
# Add compound identities for glucose and acetate
samp_peaks.loc[(samp_peaks['retention_time'] >= 15) & 
                (samp_peaks['retention_time'] <= 15.4), 'compound'] = 'phosphate'
samp_peaks.loc[(samp_peaks['retention_time'] >= 15.5) & 
               (samp_peaks['retention_time'] <= 15.9), 'compound'] = 'glucose'
samp_peaks.loc[(samp_peaks['retention_time'] >= 24.5) & 
               (samp_peaks['retention_time'] <= 30.1), 'compound'] = 'acetate'
samp_peaks.dropna(inplace=True)

#%%
# Drop unnecessary columns and save to disk
samp_peaks.drop(columns=['sample'], inplace=True)
samp_peaks.to_csv(f'{DATA_PATH}/processed/2021-03-31_NCM3722_glucose_turnover_peaks.csv', 
                 index=False)
# %%
