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


def standardize_genes(genes, name=True, COG=True, progress=True):
    """
    Given a list of genes, standardizee their names, assign their COG designation,
    and define their sector identity

    Parameters
    ----------
    genes: str, list, or 1d numpy.ndarray
        Name or array of names you wish to convert as strings. Case insensitive
    name: bool
        If True, the EcoCyc "common name" will be returned for the gene.
    COG: bool
        If True, the COG letter and class assignment will be returned for the gene.
    progress: bool
        If True, a status bar will print summarizing the status.
    
    Returns
    -------
    [names, [cog letters, cog class], sectors] : list of lists
        A standardardized list of lists of the standardized names, 
        cogs, and sectors in that order. Any category not requested (e.g. names=False),
        an empty list will be returned in its place

    """
    if __GENES__ is None:
        _load_colidict()

    # Determine the type of the supplied gene names 
    # if (type(genes) != list) | (type(genes) != np.ndarray) | (type(genes) != str):
        # raise TypeError(f'Provided argument is of type {type(genes)}! must be str, list, or np.ndarray')
    if type(genes) == str:
        genes = [genes]
    
    # Determine the iterator
    if progress:
        iterator = tqdm.tqdm(genes, desc='Converting gene names...')
    else:
        iterator = genes
    # Do the standardization
    std = [[], [[], []]]


    for g in iterator:
        try:
            sel = __GENES__[g.lower()]
            if name:
                std[0].append(sel['common_name'])
            if  COG:
                std[1][0].append(sel['cog_letter'])
                std[1][1].append(sel['cog_class'])
            # if sector:
                # std[2].append()
        except KeyError:
           warnings.warn(f'Provided gene {g} not present in master gene list.')
    else: 
        return std
        
def assign_COG_classification(genes, progress=True):
    """
    Given a gene name (or list of names) return the assigned COG category.

    Parameters
    ----------
    genes: str, list, or 1d numpy.ndarray
        Gene name or array of gene names you wish to assign cog classes. 
        Case insensitive
    progress: bool
        If True, a status bar will print summarizing the status.
    
    Returns
    -------
    genes_std: str, list
       Standardized list of gene names following the Ecocyc  

    """
    if __GENES__ is None:
        _load_colidict()

    # Determine the type of the supplied gene names 
    if type(genes) == 'str':
        genes = []
    elif (type(genes) != list) | (type(genes) != np.ndarray):
        raise TypeError(f'Provided argument is of type {type(genes)}! must be str, list, or np.ndarray')
    
    # Determine the iterator
    if progress:
        iterator = tqdm.tqdm(genes, desc='Converting gene names...')
    else:
        iterator = genes
    # Do the standardization
    std = []
    for g in iterator:
        try:
            std.append(genes[g.lower()][''])
        except KeyError:
            raise UserWarning(f'Provided gene {g} not present in master gene list.')
    # Return the correct thing 
    if len(std)  == 1:
        return std[0]
    else: 
        return std

def numeric_formatter(values, digits=3, sci=False, unit=''):
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