#%%
import numpy as np 
import pandas as pd 
import diaux.model
import diaux.viz
import imp 
imp.reload(diaux.model)
import altair as alt
from altair_saver import save
colors, palette = diaux.viz.altair_style()
alt.data_transformers.disable_max_rows()


# Deifne the constants
OD_CONV = 6E17
AVO = 6.022E23
VOL = 1E-3

# Define the quantitative parameters
gamma_max =  17.1 * 3600 / 7459
nu_max = 2.5
omega = 0.37
precursor_mass_ref = 2E-3
Km = 5E-6

def optimal_allocation(phi_X, nu, gamma, Kd, phi_O=0.5):
    numer = (phi_O + phi_X - 1) * (2 * Kd * nu * gamma - nu**2  -nu * gamma + np.sqrt(Kd * nu * gamma) * (nu - gamma))
    denom = 4 * Kd * nu * gamma - nu**2 - 2 * nu * gamma - gamma**2
    return -numer / denom

# Define allocation parameters
phi_O = 0.5
phi_X = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07])
num_muts = len(phi_X)
# Determine the optimal allocation o
phi_R = optimal_allocation(phi_X, nu_max, gamma_max, precursor_mass_ref, phi_O)
phi_P = 1 - phi_O - phi_X - phi_R

# Set up initial conditions
od_init = 0.04 / 5
M = np.ones(num_muts) * od_init * OD_CONV
Mr = phi_R * M 
Mp = phi_P * M
precursors = 4.5E-4 * M 
nutrient_conc = 0.006 #in M
nutrients = np.array([nutrient_conc * AVO * VOL])

# Define the timesteps and dilution

nutrient_mass_init = nutrients
time = np.linspace(0, 3, 200)
#
# Define the parameters and arguments
params = [M, Mr, Mp, precursors, nutrients]
flat_params = [value for param in params for value in param]

args = (gamma_max, 
            nu_max, 
            precursor_mass_ref,
            Km, 
            omega,
            phi_R, 
            phi_P, 
            num_muts)

#%%
# Specify the nutrient resets
nutrient_dict = {-1: nutrients}
colnames_anchor = ['protein_mass', 'ribosome_mass', 'metabolic_mass', 'precursors']
colnames = [f'{k}__{i}' for i in range(num_muts) for k in colnames_anchor]
colnames.append('nutrients')

# Specify the dilution parameters
dilutions = 25 
target_mass = 0.04 * OD_CONV
dilution_df = diaux.model.dilution_cycle(time, diaux.model.single_nutrient, 
                                    flat_params, args, num_dilutions=dilutions, 
                                    target_mass=target_mass,
                                    nutrient_dict=nutrient_dict,
                                    colnames=colnames,
                                    num_muts=num_muts)

#%%
shared_dfs = pd.DataFrame(['time'])
anchors = []
for col in colnames:
        _dilution_df = pd.DataFrame([])
        if '__' in col:
                name, mut = col.split('__')
                if name not in anchors:
                        anchors.append(name)
# N_muts = 5 
_dfs = []
for n in range(num_muts):
        _df = pd.DataFrame([])
        for i, name in enumerate(anchors):
                _df['time'] = dilution_df['time']
                _df['idx'] = n
                _df[name] = dilution_df[f'{name}__{n}']
        _dfs.append(_df)

unique  = pd.concat(_dfs)               
for i, x in enumerate(phi_X):
        unique.loc[unique['idx'] == i, 'phi_X'] = x


#                 _dilution_df[name] = dilution_df[col] 
#                 _dilution_df['idx'] = mut
#                 _dilution_df['time'] = dilution_df['time']
#                 shared_dfs = shared_dfs.merge(_dilution_df, on='time') 
#         else:
#                 _dilution_df[col] = dilution_df[col]
#                 _dilution_df['time'] = dilution_df['time']
#                 tidy_dfs[1].append(_dilution_df) 


# unique = pd.concat(tidy_dfs[0], join='outer', sort=False)
# unique = unique.pivot(index='idx', columns=['value'])
# # shared = pd.concat(tidy_dfs[1],sort=False)
# unique


#%%
unique['od'] = unique['protein_mass'].values / OD_CONV
base = alt.Chart(unique, width=800, height=300)
muts = base.mark_line().encode(
                x=alt.X('time:Q', title='time [hr]'),
                y=alt.Y('od:Q', title='optical density'),
                color=alt.Color('phi_X:O', title='φx',
                                scale=alt.Scale(scheme='viridis')))
total = base.mark_line(opacity=0.2, size=4, color='black').encode(
                x=alt.X('time:Q', title='time [hr]'),
                y=alt.Y('sum(od):Q', title='optical density')
)

save(muts, './steady_state_periodic_dilution.pdf')

# %%
