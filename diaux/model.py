import numpy as np 
import pandas as pd 
import scipy.integrate

def single_nutrient(params, time, gamma_max, nu_max, precursor_mass_ref, Km, 
                    omega, phi_R, phi_P, dilution=False,nutrient_mass_init=0, OD_crit=1.0,
                    dilution_factor=1E3, volume=1E-3):
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
    dilution_interval : float or None
        The time interval at which the system is diluted for another growth 
        cycle. If `None`, no dilution occurs.
    dilution_factor : positive float 
        The factor by which the system is diluted at the specified time interval. 
        This is only used if `dilution_interval` is not `None`.
    nutrient_conc_init: positive float
        The nutrient mass of the new growth medium upon dilution. This is only 
        used if `dilution_interval` is not `None`.

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

    # Determine if the dilution should occur
    if dilution:     
        if M >= (OD_CONV * OD_crit):
            M = M / dilution_factor
            Mr = Mr / dilution_factor
            Mp = Mp / dilution_factor
            precursors = precursors / dilution_factor
            nutrients = nutrient_mass_init

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


def integrate_temporal_dilutions(time, dilution_interval, dilution_factor, 
                                 fun, params, args):
    """
    Integrates a desired function with periodic dilutions and returns a
    dataframe of the complete integration. 

    NOTE: This function assumes the final 
    """
    return False