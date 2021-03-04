#%%
import numpy as np 
import scipy.integrate
import pandas as pd 
import diaux.viz 
import altair as alt 
colors, palette = diaux.viz.altair_style()


# Define the equations for growth on a single nutrient
def single_nutrient_dynamics(params, time, gamma_max, nu_max, phi_R, 
                             phi_P, theta_0, Km, Y):
    """
    Defines the complete dynamical model for growth on a single carbon source. 

    Parameters
    ----------
    params: list, of order [M, theta_a, c, M_P]
        List of the model parameters which
    """
    # Unpack the initial conditions 
    M, theta_a, c, M_P = params

    # Compute necessary quantities
    gamma = gamma_max * theta_a / (theta_a + theta_0) # Translational capacity
    nu = nu_max * c / (c + Km) # Nutritional capacity

    # Define the system of equations
    dM_dt = phi_R * M * gamma
    dMp_dt = phi_P * dM_dt
    dtheta_dt = nu * M_P - dM_dt 
    dc_dt = - nu * M_P / Y 
    return [dM_dt, dtheta_dt, dc_dt, dMp_dt]


# Define the intial conditions
gamma_max = 3600 * 17.1 / 7459 # Given values for E 
nu_max = 2.5
phi_Q = 0.35
phi_R = np.linspace(0.01, 0.25, 10)
phi_P = 1 - phi_Q - phi_R
theta_0 = 1E-3 * 20
Km = 5E-3
Y = 0.3
c = 0.010 # in mM
M = 0.0005 
M_P = phi_P * M
theta_a = 0.1 * M

# Define the time steps.
n_steps = 500
t_end = 20 
t_start = 0
dt = (t_end - t_start)  / n_steps
time = np.linspace(t_start, t_end, n_steps)

# ITerate through each phi_R and perform the integration
dfs = []
for i, phi in enumerate(phi_R):
    out = scipy.integrate.odeint(single_nutrient_dynamics, [M, theta_a, c, M_P[i]],
                                 time, args=(gamma_max, nu_max, phi_R[i], phi_P[i], 
                                             theta_0, Km, Y))
    _df = pd.DataFrame(out, columns=['dM_dt', 'dtheta_dt', 'dc_dt', 'dMp_dt'])
    _df['time'] = time
    _df['relative_dM_dt'] = _df['dM_dt'] / M
    _df['relative_dc_dt'] = _df['dc_dt'] / c
    _df['relative_dtheta_dt'] = _df['dtheta_dt'] / theta_a
    _df['phi_R'] = phi
    dfs.append(_df)

single_nutrient_df = pd.concat(dfs, sort=False)
single_nutrient_df.head()

base = alt.Chart(single_nutrient_df, width=300, height=250).encode(
        x=alt.X('time:Q', title='time [hr]'),
        color=alt.Color('phi_R:Q',  scale=alt.Scale(scheme='lighttealblue'), title='φR')
)

dM_dt = base.mark_line(size=2).encode(
        y=alt.Y('relative_dM_dt:Q', title='relative biomass, M(t)/M0')

)
dc_dt = base.mark_line(size=2).encode(
            y=alt.Y('relative_dc_dt:Q', axis=alt.Axis(format='%'), 
                    title='nutrient composition of medium')

)

dtheta_dt = base.mark_line(size=2).encode(
            y=alt.Y('relative_dtheta_dt:Q')

)
dM_dt | dc_dt | dtheta_dt
# %%
