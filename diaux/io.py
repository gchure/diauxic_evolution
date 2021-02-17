import numpy as np 
import pkg_resources 
import pandas as pd
import warnings
from io import StringIO
import pickle
import tqdm

__GENES__ = None
def _load_colidict():
    global __GENES__
    f = pkg_resources.resource_filename('diaux', 'package_data/coli_gene_dict.pkl')
    with open(f, 'rb') as file:
        __GENES__ = pickle.load(file)

def standardize_genes(genes, strain='MG1655', progress=True):
    """
    Given a list of genes, standardizee their names, assign their COG designation,
    and define their sector identity

    Parameters
    ----------
    genes: str, list, or 1d numpy.ndarray
        Name or array of names you wish to convert as strings. Case insensitive
    strain: str, either 'MG1655' or 'REL606'
        The strain of E. coli whose gene names you wish to annotate.    

    Returns
    -------
    [genes, [cog_class, cog_letter, cog_desc], sector]: list of lists
        A standardardized list of lists of the standardized names, 
        cogs, and sectors in that order.
    """
    global __GENES__
    if __GENES__ is None:
        _load_colidict()

    # Determine if the desired string is correct.
    try:
        strain = __GENES__[strain.lower()]
    except:
        raise ValueError(f"Strain '{strain}' not in database. Must be either 'MG1655' or 'REL606'")

    if type(genes) == str:
        genes = [genes]
    
    # Determine the iterator
    if progress:
        iterator = tqdm.tqdm(genes, desc='Converting gene names...')
    else:
        iterator = genes

    # Do the standardization
    std = [[], [[], [], []], []]
    for g in iterator:
        try:
            sel = strain[g.lower()]
            std[0].append(sel['common_name'])
            std[1][0].append(sel['cog_letter'])
            std[1][1].append(sel['cog_class'])
            std[1][2].append(sel['cog_desc'])
            std[2].append(sel['sector'])
            
        except KeyError:
           print(f'Warning: provided gene {g} not present in master gene list.')
           std[0].append(g)
           std[1][0].append('Not Found')
           std[1][1].append('Not Found')
           std[1][2].append('Not Found')
           std[2].append('Not Found')
    __GENES__ = None
    return {'names':std[0], 
            'cog_letters':std[1][0],
            'cog_classes':std[1][1],
            'cog_descs':std[1][2],
            'sectors':std[2]}
        
def numeric_formatter(values, digits=3, sci=True, unit=''):
    """
    Formats numbers to human-readable formats using single-letter abbreviations
    for orders of magntiude. 

    Parameters
    ----------
    values : list or nd-array of numeric values
        Value which you wish for format in a readable manner
    digits : int
        Number of significant figures to report

    Returns
    -------
    str_vals : list 
        List of formatted numbers of the same length as values.
    """
    base_powers = np.floor(np.log10(values))
    str_vals = []
    base_dict = {'p':[-15, -12], 'n':[-12, -9], 'µ':[-9, -6], 'm':[-6, -3],
                 '':[-3, 3], 'k':[3,6], 'M':[6, 9], 'B':[9, 12],
                 'T':[12,15]}
    for i, v in enumerate(values):
        base = base_powers[i]
        for _v, _k in base_dict.items():
            if (base >=_k[0]) & (base < _k[1]):
                n = _k[0]
                l = _v
        val = str(np.round(v*10**-n, decimals=digits)) 
        if len(val.replace('.', '')) <= digits:
            val = val
        else:
            if '.' in val:
                end = digits + 1
            if (val[-1] == '.'):
                end = -1
            # if val[digits-1] == '.':
                # end = digits + 1
            # if val[1] == '.':
                # end = digits + 1
            else:
                end = digits  
            val = val[:end]
        if (sci == True) & (l == 'B'):
            l = 'G'
        str_vals.append(f'{val} {l}{unit}')
    return str_vals 


def parse_HPLC_output(file_path, peak_table=True, chromatogram=True):
    """
    Parses the textfile output of the Shimadzu HPLC. Note that this function 
    assumes the desired detector output is "Detector B", that there is a "Compound"
    entry, and that the LC Chromatogram is the last entry in the file.

    TODO: Make this less hardcoded

    Parameters
    ----------
    file_path: str
        Path to the text file output
    peak_table: bool
        If True, the peak table will be identified and returned
    chromatogram: bool 
        If True, the chromatogram will be identified and returned

    Returns
    -------
    peak_table_df, chromatogram_df: pandas DataFrame  
        Pandas DataFrames for the peak table and chromatogram. Returns follow 
        the kwargs

    Raises
    ------
    TypeError:
        TypeError is returned if the file_path is not a string
    RuntimeError:
        RuntimeError is returned if neither `peak_table` nor `chromatogram` is True
    """

    # Ensure types are correct
    if type(file_path) is not str:
        raise TypeError(f'file_path must be a string. Provided argument is of type {type(file_path)}')
    if (peak_table != True) & (chromatogram != True):
        raise RuntimeError(f'No objects set to return. Both `peak_table` and `chromatogram` are False')

    # Load the file and read the lines
    with open(file_path, 'r') as f:
        lines = f.read()

    # Determine if the peak table should be isolated
    out = []
    if peak_table:
        peaks = '\n'.join(lines.split('[Peak Table(Detector B)]')[1].split('[Compound Results')[0][1:].split('\n')[1:])
        peak_table_df = pd.read_csv(StringIO(peaks))
         
        # Clean and properly annotate the peak table.
        peak_table_df = peak_table_df[['Peak#', 'R.Time', 'I.Time', 'F.Time', 'Area', 'Height']]
        peak_table_df.rename(columns={'Peak#':'peak_id', 'R.Time':'retention_time', 
                       'I.Time':'initial_time', 'F.Time':'final_time',
                       'Area':'peak_area', 'Height':'peak_height'}, inplace=True)
        out.append(peak_table_df)
    # Determine if the chromatogram should be isolated
    if chromatogram:
        # Isolate the chromatogram
        chrom =  '\n'.join(lines.split('[LC Chromatogram(Detector B-Ch1)]')[1].split('Multiplier')[1].split('\n')[1:])
        chromatogram_df = pd.read_csv(StringIO(chrom))
        chromatogram_df.rename(columns={'R.Time (min)':'time_min', 
                              'Intensity':'intensity_mV'},
                                  inplace=True)
        out.append(chromatogram_df)

    # Determine the return 
    if len(out) == 1:
        return out[0] 
    else:
        return out