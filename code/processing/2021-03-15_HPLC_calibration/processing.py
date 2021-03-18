#%%
import numpy as np
import pandas as pd 
import diaux.io
import tqdm
from io import StringIO
import glob

# Define constants
DATA_PATH = '../../../data'

# Load the list of turnover files 
cal_samples = glob.glob(f'{DATA_PATH}/2021-03-15_HPLC_calibration/raw/*.txt')

#%%    
# Generate the files for the calibration
# Instantiate storage lists for the peak tables and chromatogram
outputs = [[], []]

# Parse identifing information about the sample
for s in tqdm.tqdm(cal_samples, desc='Processing calibration measurements'):
    fname = s.split('/')[-1]
    date, _, conc = fname.split('_')    
    conc = float(conc.split('mM')[0])

    # Isolate the peak table and chromatogram
    out = diaux.io.parse_HPLC_output(s, detector='B-Ch1')

    # Add identifying information to each
    for i, df in enumerate(out):
        df['date'] = date
        df['solute_concentration_mM'] = conc
        outputs[i].append(df)

# Concatenate the dataframes
peaks = pd.concat(outputs[0])
chroms = pd.concat(outputs[1])

# Save the files to disk
peaks.to_csv(f'{DATA_PATH}/2021-03-15_HPLC_calibration/processed/2021-03-15_calibration_peaks.csv', index=False)
chroms.to_csv(f'{DATA_PATH}/2021-03-15_HPLC_calibration/processed/2021-03-15_calibration_chromatograms.csv', index=False)
