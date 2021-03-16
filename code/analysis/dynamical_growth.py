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
phi_R = 0.25 
phi_P = 0.15

# Set up initial conditions
od_init = 0.04
M = od_init * OD_CONV
Mr = phi_R * M 
Mp = phi_P * M
precursors = 4.5E-4 * M 
nutrient_conc = 0.01 #in M
nutrients = nutrient_conc * AVO * VOL

# Define the timesteps and dilution

nutrient_mass_init = nutrients
t_start = 0
t_end = 30
n_steps = 500
dt = (t_end - t_start) / n_steps
time = np.linspace(0, 30, 500)
# dilution_times = [time[100], time[200], time[400]]
# dilution_interval = 4
params = [M, Mr, Mp, precursors, nutrients]
args = (gamma_max, nu_max, precursor_mass_ref,
        Km, omega, phi_R, phi_P, False,
        nutrient_mass_init, dt)
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
