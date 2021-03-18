import numpy as np 
import pandas as pd 
import scipy.integrate
import tqdm

def single_nutrient(params, time, gamma_max, nu_max, precursor_mass_ref, Km, 
                    omega, phi_R, phi_P, volume=1E-3):
    """
    Defines the system of ordinary differenetial equations (ODEs) which describe 
    accumulation of biomass on a single nutrient source. 

    Parameters
    ----------

    params: list, [M, Mr, Mp, precursors, nutrients]
        A list of the parameters whose dynamics are described by the ODEs.
        M : positive float
            Total protein biomass of the system
        Mr : positive float, must be < M 
            Ribosomal protein biomass of the system
        Mp : positive float, must be < M
            Metabbolic protein biomass of the system 
        precursors : positive float
            Mass of precursors in the cell. This is normalized to 
            total protein biomass when calculating the translational 
            capacity.
        nutrients  : positive float
            Mass of nutrients in the system.
    time : float
        Evaluated time step of the system.
    gamma_max: positive float 
        The maximum translational capacity in units of inverse time.
    nu_max : positive float
        The maximum nutritional capacity in units of inverse time. 
    precursor_conc_ref : positive float 
        The dissociation constant of charged tRNA to the elongating ribosome.   
    Km : positive float
        The Monod constant for growth on the specific nutrient source. 
        This is in units of molar.
    omega: positive float
        The yield coefficient of the nutrient source in mass of amino acid 
        produced per mass of nutrient.
    phi_R : float, [0, 1]
        The fraction of the proteome occupied by ribosomal protein mass
    phi_P : float, [0, 1] 
        The fraction of the proteome occupied by metabolic protein mass
    volume: float, default 1 mL
        The volume of the system for calculation of concentrations.

    Returns
    -------
    out: list, [dM_dt, dMr_dt, dMp_dt, dprecursors_dt, dnutrients_dt]
        A list of the evaluated ODEs at the specified time step.

        dM_dt : The dynamics of the total protein biomass.
        dMr_dt : The dynamics of the ribosomal protein biomass.
        dMp_dt : the dynamics of the metabolic protein biomass.
        dprecursors_dt : The dynamics of the precursor/charged-tRNA pool.
        dnutrients_dt :  The dynamics of the nutrients in the growth medium
    """
    # Define constants 
    AVO = 6.022E23
    OD_CONV = 6E17
    #TODO: Put in data validation
        
    # Unpack the parameters
    M, Mr, Mp, precursors, nutrients = params

    # Compute the precursor mass fraction and nutrient concentration
    precursor_mass_frac = precursors / M
    nutrient_conc = nutrients / (AVO * volume)

    # Compute the two capacities
    gamma = gamma_max * precursor_mass_frac / (precursor_mass_frac + precursor_mass_ref)
    nu = nu_max * nutrient_conc / (nutrient_conc + Km)

    # ODEs for biomass accumulation
    dM_dt = gamma * Mr
    dMr_dt = phi_R * dM_dt
    dMp_dt = phi_P * dM_dt
    # ODE for precursors and nutrients
    dprecursors_dt = nu * Mp - dM_dt
    dnutrients_dt = -nu * Mp/ omega

    return [dM_dt, dMr_dt, dMp_dt, dprecursors_dt, dnutrients_dt]

def single_nutrient_diverse_pop(params, time, gamma_max, nu_max, precursor_mass_ref, Km, 
                    omega, phi_R, phi_P, num_muts=3, volume=1E-3):
    """
    Defines the system of ordinary differenetial equations (ODEs) which describe 
    accumulation of biomass on a single nutrient source. 

    Parameters
    ----------

    params: list of floats or arrays, [M, Mr, Mp, precursors, nutrients]
        A list of the parameters whose dynamics are described by the ODEs.
        M : array, positive float
            Total protein biomass of each mutant in the system
        Mr : array, positive float, must be < M 
            Ribosomal protein biomass of each mutant in the the system
        Mp : array, positive float, must be < M
            Metabolic protein biomass of each mutant in the system 
        precursors : array, positive float
            Mass of precursors for each mutant in the system. This is normalized 
            to total protein biomass when calculating the translational 
            capacity.
        nutrients  : positive float
            Mass of nutrients in the system. This is consumed by each mutant.
    time : float
        Evaluated time step of the system.
    gamma_max: array, positive float 
        The maximum translational capacity for each mutant in units of inverse 
        time.
    nu_max : array, positive float
        The maximum nutritional capacity for each mutant in units of inverse 
        time. 
    precursor_conc_ref : array, positive float 
        The dissociation constant of charged tRNA to the elongating ribosome 
        for each mutant.  
    Km : array, positive float
        The Monod constant for growth on the specific nutrient source for each 
        mutant in the system. This is in units of molar.
    omega: positive float
        The yield coefficient of the nutrient source in mass of amino acid 
        produced per mass of nutrient.
    phi_R : array, float, [0, 1]
        The fraction of the proteome occupied by ribosomal protein mass for 
        each mutant in the system. 
    phi_P : array, float, [0, 1] 
        The fraction of the proteome occupied by metabolic protein mass for 
        each mutant in the system.
    volume: float, default 1 mL
        The volume of the system for calculation of concentrations.

    Returns
    -------
    out: list, [dM_dt, dMr_dt, dMp_dt, dprecursors_dt, dnutrients_dt]
        A list of the evaluated ODEs at the specified time step.
        dM_dt : The dynamics of the total protein biomass.
        dMr_dt : The dynamics of the ribosomal protein biomass.
        dMp_dt : the dynamics of the metabolic protein biomass.
        dprecursors_dt : The dynamics of the precursor/charged-tRNA pool.
        dnutrients_dt :  The dynamics of the nutrients in the growth medium
    """
    # Define constants 
    AVO = 6.022E23
    OD_CONV = 6E17
    #TODO: Put in data validation
        
    # Unpack the parameters
    nutrients = params[-1]
    M, Mr, Mp, precursors = np.reshape(params[:-1], (4, num_muts))
    
    # Compute the precursor mass fraction and nutrient concentration
    precursor_mass_frac = precursors / M
    nutrient_conc = nutrients / (AVO * volume)

    # Compute the two capacities
    gamma = gamma_max * precursor_mass_frac / (precursor_mass_frac + precursor_mass_ref)
    nu = nu_max * nutrient_conc / (nutrient_conc + Km)

    # ODEs for biomass accumulation
    dM_dt = gamma * Mr
    dMr_dt = phi_R * dM_dt
    dMp_dt = phi_P * dM_dt

    # ODE for precursors and nutrients
    dprecursors_dt = nu * Mp - dM_dt
    dnutrients_dt = -np.sum(nu * Mp/ omega)

    # Reshape for the output
    _out = [dM_dt, dMr_dt, dMp_dt, dprecursors_dt]
    out = [value for deriv in _out for value in deriv]
    out.append(dnutrients_dt)
    return out

def temporal_dilutions(time, fun, fun_params, fun_args, nutrient_dict, num_dilutions=10, 
                        dilution_factor=1E3, colnames=None, **int_kwargs):
    """
    Integrates a desired function with periodic dilutions and returns a
    dataframe of the complete integration. 

    Parameters 
    -----------
    time: numpy-array
        The time interval to integrate for a single growth cycle. This 
        time interval will be repeated for each dilution. 
    fun: function
        The function you wish to integrate
    fun_params : list
        List of parameters to feed into the function
    fun_args : tuple 
        Arguments to feed the integration function
    nutrient_dict : dict
        A dictionary of the indices and values to reset the nutrient conditions 
        for each dilution event. The keys correspond to the indices those of the 
        `fun_params` which define the nutrient conditions. The value corresponds 
        to the desired reset value. 
    num_dilutions : int 
        The number of dilution cycles that should be performed
    dilution_factor : float or int
        The factor by which the parameters should be decreased for each dilution 
        event. Note that this does not apply to the nutrient parameters which 
        are reset by `nutrient_dict.`
    colnames : list of str, optional
        The desired column names of the output. If `None`, columns will be 
        left arbitrarily named.
    **int_kwargs: dict
        kwargs to be fed to the ODE solver.
    """

    # TODO: Put in type checks.

    # Perform the initial integration
    out = scipy.integrate.odeint(fun, fun_params, time, args=fun_args,
                                         **int_kwargs)
    # Insstantiate the dataframes
    if colnames != None:
        initial_df  = pd.DataFrame(out, columns=colnames)
    else:
        initial_df = pd.DataFrame(out)

    # Specify that this is the inaugural integration
    initial_df['dilution_cycle'] = 0

    # Keep track of the time 
    initial_df['time'] = time

    stop_time = time[-1]
    
    # Add to the storage list
    dfs = [initial_df]

    # Iterate through each dilution cycle.
    for n in tqdm.tqdm(range(num_dilutions)):
        # Reset the parameters
        fun_params = list(out[-1, :] / dilution_factor)

        for k, v in nutrient_dict.items():
            fun_params[k] = v
        
        # Rerun the integration
        out = scipy.integrate.odeint(fun, fun_params, time, args=fun_args,
                                         **int_kwargs)

        # Set up the dataframe 
        if colnames != None:
            _df = pd.DataFrame(out, columns=colnames)
        else:
            _df = pd.DataFrame(out)

        # Update time and keep track of the dilution cycle
        _df['time'] = stop_time + time
        _df['dilution_cycle'] = n + 1

        # Add the dataframe to the storage list
        dfs.append(_df)
        
        # Update the stoptime
        stop_time += time[-1]

    # Concatenate the dataframes and return
    dil_df = pd.concat(dfs)
    return  dil_df