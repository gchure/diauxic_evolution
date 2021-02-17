#%%
import numpy as np
import pandas as pd 
import diaux.io
import tqdm
from io import StringIO
import glob

# Load the list of turnover files 
growth_samples = glob.glob('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/raw/*HPLC*diauxic*.txt')
cal_samples = glob.glob('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/raw/*HPLC*calibration*.txt')


# Instantiate storage lists for the peak tables and chromatogram
outputs = [[], []]

# Parse identifing information about the sample
for s in tqdm.tqdm(growth_samples, desc='Processing diauxic shift measurements'):
    fname = s.split('/')[-1]
    date, _, strain, _, _, _, sample_id = fname.split('_')
    sample_id = int(sample_id.split('.')[0])
    sample_id

    if strain == 'NCM':
        strain = 'NCM3722'

    # Isolate the peak table and chromatogram
    out = diaux.io.parse_HPLC_output(s)

    # Add identifying information to each
    for i, df in enumerate(out):
        df['strain'] = strain
        df['date'] = date
        df['sample_id'] = sample_id
        outputs[i].append(df)

# Concatenate the dataframes
peaks = pd.concat(outputs[0])
chroms = pd.concat(outputs[1])

# Save the files to disk
peaks.to_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_diauxic_shift_peaks.csv', index=False)
chroms.to_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_diauxic_shift_chromatograms.csv', index=False)

#%%
# Generate the files for the calibration
# Instantiate storage lists for the peak tables and chromatogram
outputs = [[], []]

# Parse identifing information about the sample
for s in tqdm.tqdm(cal_samples, desc='Processing calibration measurements'):
    fname = s.split('/')[-1]
    date, _, _, carbon, conc = fname.split('_')
    conc = conc.split('mM')[0]

    # Isolate the peak table and chromatogram
    out = diaux.io.parse_HPLC_output(s)

    # Add identifying information to each
    for i, df in enumerate(out):
        df['date'] = date
        df['carbon_source'] = carbon
        df['concentration_mM'] = conc
        outputs[i].append(df)

# Concatenate the dataframes
peaks = pd.concat(outputs[0])
chroms = pd.concat(outputs[1])

# Save the files to disk
peaks.to_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_calibration_peaks.csv')
chroms.to_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_calibration_chromatograms.csv')
# %%
