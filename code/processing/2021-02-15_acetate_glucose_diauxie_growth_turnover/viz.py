#%%
import pandas as pd 
import altair as alt 
from altair_saver import save
import diaux.viz
colors, palette = diaux.viz.altair_style()

# Load the complete data data = pd.read_csv('../../../data/2021-02-15_NCM3722_glucose_acetate_diauxic_shift/processed/2021-02-15_glucose_acetate_turnover_growth.csv') 
# Set up the chart of the growth curve. 
growth = alt.Chart(data).mark_point(size=80).encode(
                x=alt.X('elapsed_time_hr:Q', title='elapsed time [hr]'),
                y=alt.Y('od_600nm:Q', title='optical density [a.u.]'),
                color=alt.Color('region:N', title='growth phase'))

# Set up the chart of the glucose concentration
_data = data[data['glucose_concentration_mM'] >= 0]
gluc = alt.Chart(_data).encode(
            x=alt.X('elapsed_time_hr:Q', title='elapsed time [hr]'),
            y=alt.Y('glucose_concentration_mM:Q', title='glucose [mM]', 
                axis=alt.Axis(titleColor=colors['primary_purple']))
)
gluc_point = gluc.mark_point(color=colors['primary_purple'], size=80)
gluc_line = gluc.mark_line(color=colors['primary_purple'], size=2)
ace = alt.Chart(_data).encode(
            x=alt.X('elapsed_time_hr:Q', title='elapsed time [hr]'),
            y=alt.Y('acetate_concentration_mM:Q', title='acetate [mM]',
                    axis=alt.Axis(titleColor=colors['primary_red']))
)
ace_point = ace.mark_point(color=colors['primary_red'], size=80)
ace_line = ace.mark_line(color=colors['primary_red'], size=2)

gluc_plot = (growth + gluc_line + gluc_point).resolve_scale(y='independent')
ace_plot = (growth + ace_line + ace_point).resolve_scale(y='independent')

save(gluc_plot, './output/2021-02-15_glucose_growth_turnover.pdf')
save(ace_plot, './output/2021-02-15_glucose_acetate_turnover.pdf')
