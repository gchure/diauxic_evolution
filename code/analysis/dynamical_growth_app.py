#%%
import numpy as np 
import pandas as pd 
import tqdm
# import bokeh.io 
# import bokeh.plotting
# import bokeh.widgets
import scipy.integrate
import diaux.viz
import altair as alt
colors, palette = diaux.viz.altair_style()
#%%
# Define the conversion factor for OD to AA unit
OD_CONV = 6E17 # in AA / mL at OD=1
AVO = 6.022E23 # Avogadros number
VOL = 1E-3 # System volume in L
# Define the functions for single and dual nutrient growth
def single_nutrient(params, t, gamma_max, nu_max, theta_0, K_m, phi_R, phi_P, 
                    omega_x, vol):
    M, Mr, Mp, m_x, n_x = params

    # Compute the nutrient mass fraction
    theta =  m_x / M
    c_x = (n_x / (AVO * vol))

    # Compute the translational capacity
    gamma = gamma_max * theta / (theta + theta_0)
    nu = nu_max * c_x / (c_x + K_m)

    # Describe biomass accumulation
    dM_dt = gamma * Mr 
    dMr_dt =  phi_R * dM_dt
    dMp_dt = phi_P * dM_dt 

    # Describe amino acid concentrations
    dmx_dt = nu * Mp - dM_dt
    dnx_dt = -nu * Mp / omega_x
    return [dM_dt, dMr_dt, dMp_dt, dmx_dt, dnx_dt]

# Constants
gamma_max = (17.1 * 3600) / 7459 # in s^-1 using numbers for E. coli
omega_x = 0.3 # about 1 AA for every three nutrient molecules
nu_max =  2.5 # (200 * omega_x * 3600) / (1E3 * 1E2)
K_m = 5E-6 # Monod constant for growth on glucose
theta_0 = 0.002

# Set up allocation parameters, assuming they don't change with concentration
phi_O = 0.5
phi_R = 0.2 
phi_P = 1 - phi_O - phi_R

# Set up initial conditions
od_init = 0.04 # in OD600nm units
M = od_init * OD_CONV# in units of AAs per unit volume
Mr = phi_R * M # Number of AAs in ribosome mass
Mp = phi_P * M # Number of AAs in metabolic mass

# Mass of charged tRNAs
m_x = 0.1 * M 

# Nutrient concentration and numbers
c_x = 0.61E-3 # Conceentration in M
n_x = c_x * AVO * VOL

# Define time step for the integration
n_steps = 500
time_start = 0
time_end = 10 
time = np.linspace(time_start, time_end, n_steps)

# Pack up parameters and arguments
_params = [M, Mr, Mp, m_x, n_x]
_args = (gamma_max, nu_max, theta_0, K_m, phi_R, phi_P,
        omega_x, VOL)

# Define the time steps
out = scipy.integrate.odeint(single_nutrient, _params, time, args=_args)

df = pd.DataFrame(out, columns=['dM_dt', 'dMr_dt', 'dMp_dt', 'dmx_dt', 'dnx_dt'])
df['time'] = time
df['theta'] = df['dmx_dt'].values / df['dM_dt'].values
df['phi_R'] = df['dMr_dt'].values / df['dM_dt'].values
df['phi_P'] = df['dMp_dt'].values / df['dM_dt'].values
df['dcx_dt'] = df['dnx_dt'].values / (AVO * VOL)
df['od'] = df['dM_dt'].values /(OD_CONV)
alt.Chart(df).mark_line().encode(x='time:Q', y=alt.Y('od:Q'))

# %%

def diauxic_nutrient(params, t, gamma_max, nu_x_max, nu_y_max, theta_0, K_mx, K_my, 
                    phi_x_max, phi_y_max, phi_R, omega_x, omega_y,  vol):
    M, Mr, Mx, My, m, n_x, n_y = params

    # Compute the nutrient mass fraction
    theta =  m / M
    c_x = (n_x / (AVO * vol))
    c_y = (n_y / (AVO * vol))

    # Compute the translational capacity
    gamma = gamma_max * theta / (theta + theta_0)
    nu_x = nu_x_max * c_x / (c_x + K_mx)
    nu_y = nu_y_max * c_y / (c_y + K_my)

    # Define the allocation     
    phi_X = phi_x_max * (c_x / (c_x + K_mx))  
    phi_Y = phi_y_max * (1 - (c_x / (c_x + K_mx)))
    # phi_R = 1 - (phi_O + phi_X + phi_Y)

    # Describe biomass accumulation
    dM_dt = gamma * Mr 
    dMr_dt = phi_R * dM_dt
    dMx_dt = phi_X * dM_dt 
    dMy_dt = phi_Y * dM_dt 

    # Describe amino acid concentrations
    dm_dt = nu_x * Mx + nu_y * My - dM_dt
    dnx_dt = -nu_x * Mx / omega_x
    dny_dt = -nu_y * My / omega_y
    return [dM_dt, dMr_dt, dMx_dt, dMy_dt, dm_dt, dnx_dt, dny_dt]

# Constants
gamma_max = (17.1 * 3600) / 7459 # in s^-1 using numbers for E. coli
omega_x = 0.3 # about 1 AA for every three nutrient molecules
omega_y = 0.1 # about 1 AA for every ten nutrient molecules
nu_x_max = 2.5 # (200 * omega_x * 3600) / (1E3 * 1E2)
nu_y_max = 5            
K_mx = 5E-6 # Monod constant for growth on glucos
K_my = 5E-6 #
theta_0 =0.002 

# Set up allocation parameters, assuming they don't change with concentration
phi_O = 0.6
phi_x = 0.2
phi_y = 0.1 
phi_R = 0.2

# Set up initial conditions
od_init = 0.04 # in OD600nm units
M = od_init * OD_CONV# in units of AAs per unit volume
Mx = phi_x * M
My = 0
Mr = M - Mx - My - phi_O * M

# Mass of charged tRNAs
m = 0.1 * theta_0 * M 

# Nutrient concentration and numbers
c_x = 0.2E-3 # Conceentration in M
c_y = 100E-3# Conceentration in M
n_x = c_x * AVO * VOL
n_y = c_y * AVO * VOL

# Define time step for the integration
n_steps = 500
time_start = 0
time_end =  10 
time = np.linspace(time_start, time_end, n_steps)

# Pack up parameters and arguments
_params = [M, Mr, Mx, My, m, n_x, n_y]
_args = (gamma_max, nu_x_max, nu_y_max, theta_0, K_mx, 
         K_my, phi_x, phi_y, phi_R, omega_x, omega_y, VOL)
# Define the time steps
out = scipy.integrate.odeint(diauxic_nutrient, _params, time, args=_args)

df = pd.DataFrame(out, columns=['dM_dt', 'dMr_dt', 'dMx_dt', 'dMy_dt', 'dm_dt', 'dnx_dt', 'dny_dt'])
df['time'] = time
df['theta'] = df['dm_dt'].values / df['dM_dt'].values
df['phi_R'] = df['dMr_dt'].values / df['dM_dt'].values
df['phi_x'] = df['dMx_dt'].values / df['dM_dt'].values
df['phi_y'] = df['dMy_dt'].values / df['dM_dt'].values
df['dcx_dt'] = df['dnx_dt'].values / (AVO * VOL)
df['dcy_dt'] = df['dny_dt'].values / (AVO * VOL)
df['od'] = df['dM_dt'].values /(OD_CONV)
base = alt.Chart(df)
od = base.mark_line().encode(x='time:Q', y='od:Q')
# phi_x = base.mark_line().encode(x='time:Q', y='phi_x:Q')
# phi_y = base.mark_line().encode(x='time:Q', y='phi_y:Q')
# od & phi_y
#%%
phi_y = np.linspace(0.01, 0.18, 20)
# Define time step for the integration
n_steps = 200
time_start = 0
time_end = 10 
time = np.linspace(time_start, time_end, n_steps)

# Pack up parameters and arguments
dfs = []
for i, phi in enumerate(tqdm.tqdm(phi_y)):
    _params = [M, Mr, Mx, My, m, n_x, n_y]
    _args = (gamma_max, nu_x_max, nu_y_max, theta_0, K_mx, 
             K_my, phi_x, phi, phi_R, omega_x, omega_y, VOL)
    # Define the time steps
    out = scipy.integrate.odeint(diauxic_nutrient, _params, time, args=_args)

    df = pd.DataFrame(out, columns=['dM_dt', 'dMr_dt', 'dMx_dt', 'dMy_dt', 'dm_dt', 'dnx_dt', 'dny_dt'])
    df['time'] = time
    df['theta'] = df['dm_dt'].values / df['dM_dt'].values
    df['phi_R'] = df['dMr_dt'].values / df['dM_dt'].values
    df['phi_x'] = df['dMx_dt'].values / df['dM_dt'].values
    df['phi_y'] = df['dMy_dt'].values / df['dM_dt'].values
    df['dcx_dt'] = df['dnx_dt'].values / (AVO * VOL)
    df['dcy_dt'] = df['dny_dt'].values / (AVO * VOL)
    df['od'] = df['dM_dt'].values /(OD_CONV)
    df['phi_y_max'] = phi
    dfs.append(df)

titrate_phi_df = pd.concat(dfs, sort=False)
#%%
from altair_saver import save
base = alt.Chart(titrate_phi_df)
od = base.mark_line(clip=True).encode(
            x=alt.X('time:Q', title='time [hr]'),
            y=alt.Y('od:Q', title='approximate optical density [a.u.]',
                    scale=alt.Scale(domain=[0.03, 0.8], type='log')),
            color=alt.Color('phi_y_max:Q', title='allocation φy'))

# # # plot the allocation   
# phi_R = base.mark_line(clip=True).encode(
#             # x=alt.X('time:Q', title='time [hr]'),
#             y=alt.Y('phi_R:Q', title='allocation φR'),
#             color=alt.Color('phi_y_max:Q', title='allocation φy'))
save(od, 'log_diauxic_optical_density.png')
od
# %%
base.mark_line(clip=True).encode(
        x=alt.X('time:Q', title='time [hr]'),
        y=alt.Y('theta:Q'),
        color=alt.Color('phi_y_max:Q')
)
# %%
