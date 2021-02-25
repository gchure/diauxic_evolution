#%%
import numpy as np
import pandas as pd 
import diaux.io
import tqdm
from io import StringIO
import glob

# Define constants
DATE = "2021-02-22"
STRAIN = 'NCM3722'
CONDITION = 'glucose_glucose+acetate_turnover'
DATA_PATH = '../../../data'
# Load the list of turnover files 
growth_samples = glob.glob(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}/raw/*HPLC*_r*_t*.txt')
cal_samples = glob.glob(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}/raw/*HPLC*_calibration*.txt')

#%%    
# Instantiate storage lists for the peak tables and chromatogram
outputs = [[], []]

# Parse identifing information about the sample
for s in tqdm.tqdm(growth_samples, desc='Processing turnover measurements'):
    fname = s.split('/')[-1]
    date, _, strain, _, glucose_conc, _, acetate_conc, replicate, timepoint = fname.split('_')
    glucose_conc = float(glucose_conc.split('mM')[0])
    acetate_conc = float(acetate_conc.split('mM')[0])
    replicate = int(replicate[1])
    timepoint = int(timepoint[1])

    # Isolate the peak table and chromatogram
    out = diaux.io.parse_HPLC_output(s)

    # Add identifying information to each
    for i, df in enumerate(out):
        df['strain'] = strain
        df['date'] = date
        df['replicate'] = replicate
        df['sample_id'] = timepoint
        df['glucose_mM'] = glucose_conc
        df['acetate_mM'] = acetate_conc
        outputs[i].append(df)

# Concatenate the dataframes
peaks = pd.concat(outputs[0])
chroms = pd.concat(outputs[1])

# Save the files to disk
peaks.to_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}/processed/{DATE}_{STRAIN}_{CONDITION}_peaks.csv', index=False)
chroms.to_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}/processed/{DATE}_{STRAIN}_{CONDITION}_chromatograms.csv', index=False)

#%%
# Generate the files for the calibration
# Instantiate storage lists for the peak tables and chromatogram
outputs = [[], []]

# Parse identifing information about the sample
for s in tqdm.tqdm(cal_samples, desc='Processing calibration measurements'):
    fname = s.split('/')[-1]
    date, _, strain, _, glucose_conc, _, acetate_conc, _ = fname.split('_')
    glucose_conc = float(glucose_conc.split('mM')[0])
    acetate_conc = float(acetate_conc.split('mM')[0])

    # Isolate the peak table and chromatogram
    out = diaux.io.parse_HPLC_output(s)

    # Add identifying information to each
    for i, df in enumerate(out):
        df['date'] = date
        df['glucose_mM'] = glucose_conc
        df['acetate_mM'] = acetate_conc
        outputs[i].append(df)

# Concatenate the dataframes
peaks = pd.concat(outputs[0])
chroms = pd.concat(outputs[1])

# Save the files to disk
peaks.to_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}/processed/{DATE}_{STRAIN}_calibration_peaks.csv', index=False)
chroms.to_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}/processed/{DATE}_{STRAIN}_calibration_chromatograms.csv', index=False)



# %%
