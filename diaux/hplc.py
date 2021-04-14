#%%
import pandas as pd

def assign_peaks(peaks, dropna=True):
    """
    Assign compound identities to peaks based off of standard retention times.
    """
    _peaks = peaks.copy(deep=True)
    _peaks.loc[(peaks['retention_time'] >=10) & (peaks['retention_time'] <= 11.5),
            'compound'] = 'NaCl'
    _peaks.loc[(peaks['retention_time'] >=13.5) & (peaks['retention_time'] <= 14.7),
            'compound'] = 'phosphate'

    _peaks.loc[(peaks['retention_time'] >=24) & (peaks['retention_time'] <= 26),
            'compound'] = 'acetate'
    _peaks.loc[(peaks['retention_time'] >=21.9) & (peaks['retention_time'] <= 22.3),
            'compound'] = 'glycerol'
    _peaks.loc[(peaks['retention_time'] >=15) & (peaks['retention_time'] <= 16),
            'compound'] = 'glucose'
    if dropna:
        _peaks.dropna(inplace=True)
    return _peaks

