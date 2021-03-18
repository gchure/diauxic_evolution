#%%
import numpy as np 
import pandas as pd 
import scipy.integrate 
import diaux.model
import diaux.viz
import imp 
imp.reload(diaux.model)
import altair as alt
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

# Define allocation parameters
phi_O = 0.5
phi_R = np.array([0.25, 0.22, 0.21, 0.15, 0.1])
phi_P = 1 - phi_O - phi_R

# Set up initial conditions
od_init = 0.04 / 5
M = np.ones(5) * od_init * OD_CONV
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

args = (gamma_max, nu_max, precursor_mass_ref,
        Km, omega, phi_R, phi_P, 5)

#%%

# Specify the nutrient resets
nutrient_dict = {-1: nutrients}
colnames_anchor = ['protein_mass', 'ribosome_mass', 'metabolic_mass', 'precursors']
colnames = [f'{k}__{i}' for i in range(5) for k in colnames_anchor]
colnames.append('nutrients')

# Specify the dilution parameters
dilutions = 10 
factor = 10
dilution_df = diaux.model.temporal_dilutions(time, diaux.model.single_nutrient_diverse_pop, 
                                    flat_params, args, num_dilutions=dilutions, 
                                    dilution_factor=factor,
                                    nutrient_dict=nutrient_dict,
                                    colnames=colnames)
#%%

shared_dfs = pd.DataFrame(['time'])
anchors = []
for col in colnames:
        _dilution_df = pd.DataFrame([])
        if '__' in col:
                name, mut = col.split('__')
                if name not in anchors:
                        anchors.append(name)
N_muts = 5 
_dfs = []
for n in range(N_muts):
        _df = pd.DataFrame([])
        for i, name in enumerate(anchors):
                _df['time'] = dilution_df['time']
                _df['idx'] = n
                _df[name] = dilution_df[f'{name}__{n}']
        _dfs.append(_df)

unique  = pd.concat(_dfs)               


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
alt.Chart(unique).mark_line().encode(
                x=alt.X('time:Q', title='time [hr]'),
                y=alt.Y('od:Q', title='optical density'),
                color=alt.Color('idx:N', title='mutant ID')
)
#%%
base = alt.Chart(unique).encode(
        x=alt.X('time:Q', title='time [hr]'))
od = base.mark_line().encode(
               y=alt.Y('od:Q', title='approximate optical density [a.u.]'),
               color=alt.Color('idx:N', title='mutant idx'))
od
# dilution_df['precursor_conc'] = dilution_df['precursors'].values / dilution_df['protein_mass']
# dilution_df['nutrient_conc'] = dilution_df['nutrients'].values / (AVO * VOL)
# base = alt.Chart(dilution_df).encode(
#         x=alt.X('time:Q', title='time [hr]'))
# od = base.mark_line().encode(y=alt.Y('od:Q', title='approximate optical density [a.u.]'))
# theta = base.mark_line().encode(y=alt.Y('precursor_conc:Q', title='precursor mass fraction'))
# nuts = base.mark_line().encode(y=alt.Y('nutrient_conc:Q', title='nutrient concentration [M]'))
# od & theta & nuts


#%%

out = scipy.integrate.odeint(diaux.model.single_nutrient,
                             params, time, args=args)

#%%
df = pd.DataFrame(out, columns=['biomass', 'ribosomal_mass', 'metabolic_mass', 
                                'precursors', 'nutrients'])
df['time'] = time
df['od'] = df['biomass'] / OD_CONV
df['precursor_mass_frac'] = df['precursors'].values / df['biomass']
df['nutrient_concentration'] = df['nutrients'].values / (AVO * VOL)

base = alt.Chart(df).encode(
        x=alt.X('time:Q', title='time [hr]'))
od = base.mark_line().encode(y=alt.Y('od:Q', title='approximate optical density [a.u.]'))
theta = base.mark_line().encode(y=alt.Y('precursor_mass_frac:Q', title='precursor mass fraction'))
nuts = base.mark_line().encode(y=alt.Y('nutrient_concentration:Q', title='nutrient_concentration'))
od & theta & nuts

# %%

# %%
