#%%
import numpy as np 
import scipy.integrate
import pandas as pd 
import diaux.viz 
import altair as alt 
from altair_saver import save
colors, palette = diaux.viz.altair_style()
alt.data_transformers.disable_max_rows()

# Define the equations for growth on a single nutrient
def dual_nutrient_dynamics(params, time, gamma_max, nu_x_max, nu_y_max, phi_R, phi_x_max,
                             phi_y_max, phi_y_steady, theta_0, Km_x, Km_y, omega_x, 
                             omega_y):
    """
    Defines the complete dynamical model for growth on a single carbon source. 

    Parameters
    ----------
    params: list, of order [M, theta_a, c]
        List of the model parameters which
    """
    # Unpack the initial conditions 
    M, M_x, M_y, theta_a, c_x, c_y = params  

    # Compute necessary quantities
    gamma = gamma_max * theta_a / (theta_a + theta_0) # Translational capacity
    nu_x = nu_x_max * c_x / (c_x + Km_x) # Nutritional capacity
    nu_y = nu_y_max * c_y / (c_y + Km_y) # Nutritional capacity

    # Define dependencies 
    # ALlocation for glucose
    phi_x = phi_x_max * c_x / (c_x + Km_x)

    # M_y = (1 - phi_O - phi_x - phi_R) * M
    phi_y = phi_y_max * (1 - (1 + Km_x/c_x)**-1) #* (1 - (1 + phi_y_steady/(M_y/M))**-1)


    # Define the system of equations
    dM_dt = phi_R * M * gamma

    dMx_dt = phi_x * dM_dt
    dMy_dt = phi_y * dM_dt

    # dMp_dt = (phi_R_max - phi_R) * dM_dt
    dtheta_dt = M * (nu_x * phi_x + nu_y * phi_y) - dM_dt 
    dcx_dt = -nu_x * phi_x * M / omega_x
    dcy_dt = -nu_y * phi_y * M / omega_y
    return [dM_dt, dMx_dt, dMy_dt, dtheta_dt, dcx_dt, dcy_dt] 


omega_x = 0.3
omega_y = 0.2
Km_x = 5E-6 
Km_y = 5E-6 
c_x = 0.005
c_y = 0.05
phi_R = 0.2

gamma_max = (17.1 * 3600) / 7459
nu_x_max = 2.5
nu_y_max = 2
phi_x_max = 0.3
phi_y_max = 0.02
phi_y_steady = 0.01
M = 0.001
theta_a = 0.1 * M
theta_0 = 0.001 * 20
n_steps = 500
t_start = 0
t_end = 50
dt = (t_start - t_end) / n_steps 
time = np.linspace(t_start, t_end, n_steps)
out = scipy.integrate.odeint(dual_nutrient_dynamics, [M, theta_a, c_x, c_y],
                             time, args=(gamma_max, nu_x_max, nu_y_max, phi_R,
                                         phi_x_max, phi_y_max, phi_y_steady,theta_0, Km_x, Km_y, 
                                         omega_x, omega_y))


df = pd.DataFrame(out, columns=['dM_dt', 'dtheta_dt', 'dcx_dt', 'dcy_dt'])
df['time'] = time
df['rel_m'] = df['dM_dt'].values / M
df['rel_m'] = df['rel_m'].values.astype(float)
df['rel_cx'] = df['dcx_dt'].values / c_x
df['rel_cy'] = df['dcy_dt'].values / c_y
df['theta_cell'] = df['dtheta_dt'].values / df['dM_dt'].values

base = alt.Chart(df).encode(x='time:Q')

dm = base.mark_line().encode(y='rel_m:Q')
dcx = base.mark_line().encode(y='rel_cx:Q')
dcy = base.mark_line(color=colors['primary_red']).encode(y='rel_cy:Q')
dmx = base.mark_line().encode(y='dMy_dt:Q')
dtheta = base.mark_line().encode(y='theta_cell:Q')
dm | (dcx + dcy) & dtheta

#%%

# dmx
# %%
dcy.interactive()

# %%
