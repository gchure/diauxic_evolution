import numpy as np 
import pkg_resources 
import warnings
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
