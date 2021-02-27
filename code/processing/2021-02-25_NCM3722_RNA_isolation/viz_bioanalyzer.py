#%%
import numpy as np 
import pandas as pd 
import diaux.viz
import altair as alt
from altair_saver import save
import scipy.signal
import scipy.interpolate
colors, _ = diaux.viz.altair_style()

# Define  constants of the experiment 
DATE = '2021-02-25'
STRAIN = 'NCM3722'
CONDITION = 'glucose'
OTHER = 'RNA_isolation'
DATA_PATH = '../../../data'
RIN = 10 # Annotated RIN score

data = pd.read_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}_{OTHER}/raw/{DATE}_Bioanalyzer_{STRAIN}_total_RNA.csv', comment='#')
ladder = pd.read_csv(f'{DATA_PATH}/{DATE}_{STRAIN}_{CONDITION}_{OTHER}/raw/{DATE}_Bioanalyzer_ladder.csv', comment='#')

# Identify the peaks  in the ladder
ladder_order = [25, 200, 500, 1000, 2000, 4000]


# Generate the simple plot
locs = scipy.signal.find_peaks(ladder['Value'], height=8)
peaks = ladder[ladder['Value'].isin(locs[1]['peak_heights'])]

# Add the peak identifier
peaks['size_nt'] = ladder_order
peaks['label'] = [f'{l} nt' for l in ladder_order]

# Do cubic interpolation for the standard curve
spl = scipy.interpolate.CubicSpline(peaks['Time'].values, np.log(peaks['size_nt'].values),
                                        bc_type='natural')
time_range = np.linspace(20, 60, 200)
size_nt = spl(time_range)

# Create the interpolation dataframe
interp_df = pd.DataFrame([])
interp_df['time_s'] = time_range
interp_df['size_nt'] = np.exp(size_nt)

# Generate a plot of the interpolated standard curve
points = alt.Chart(peaks).mark_point(size=100, color=colors['primary_red']).encode(
        x=alt.X('Time:Q', title='migration time [s]', scale=alt.Scale(domain=[17, 70])),
        y=alt.Y('size_nt:Q', title='fragment size [nt]', scale=alt.Scale(type='log')))
interp = alt.Chart(interp_df).mark_line(size=1, color=colors['primary_red']).encode(
        x=alt.X('time_s:Q', title='migration time [s]'),
        y=alt.Y('size_nt:Q', title='fragment size [nt]')
) 
std_curve = (points + interp).properties(title='calibration curve with cubic spline')

# Generate a plot the ladder
ladder_curve = alt.Chart(ladder).mark_line(color=colors['primary_red']).encode(
                x=alt.X('Time:Q', title='migration time [s]'),
                y=alt.Y('Value:Q', title='fluorescence [a.u.]')
)
peak_points = alt.Chart(peaks).mark_point(size=100, color=colors['primary_red']).encode(
                x=alt.X('Time:Q', title='migration time [s]'),
                y=alt.Y('Value:Q', title='fluorescence [a.u.]',
                        scale=alt.Scale(domain=[0, 30]))
)

peak_label = peak_points.mark_text(baseline='middle', dy=-10).encode(
                text='label'
)
ladder_plot = (ladder_curve + peak_points + peak_label).properties(title='ladder electropherogram')

# Assemble and save the calibration plot
cal_plot = ladder_plot & std_curve

for ext in ['pdf', 'png']:
        save(cal_plot, f'./output/{DATE}_Bioanalyzer_calibration.{ext}')

#%%

# Use the cubic interpolation to add in the approximate size for the actual 
# analysis
data.rename(columns={'Time':'migration_time_s', 'Value':'fluorescence_au'}, 
            inplace=True)
data['size_nt'] = np.exp(spl(data['migration_time_s'].values))
data = data[data['migration_time_s'] < 55]

# Identify the peaks of my sample
locs = scipy.signal.find_peaks(data['fluorescence_au'], height=18)
peaks = data[data['fluorescence_au'].isin(locs[1]['peak_heights'])]

# Add the peak identifier
peaks['label'] = [f'{int(p)} nt' for p in peaks['size_nt'].values]


# Make the plot 
sample_plot = alt.Chart(data).mark_line(clip=True, color=colors['primary_red']).encode(
                x=alt.X('size_nt:Q', title='fragment size [nt]', 
                        scale=alt.Scale(type='log', domain=[10, 4000])),
                y=alt.Y('fluorescence_au:Q', title='fluorescence [a.u.]',
                       scale=alt.Scale(domain=[0, 135]))
)

sample_peaks = alt.Chart(peaks).mark_point(size=100, color=colors['primary_red']).encode(
                x=alt.X('size_nt:Q', title='fragment size [nt]'),
                y=alt.Y('fluorescence_au:Q', title='fluorescence [a.u.]')
)

peak_labels = sample_peaks.mark_text(baseline='middle', dy=-10).encode(
                text='label'
)

sample = (sample_plot + sample_peaks + peak_labels).properties(title=f'{STRAIN} total RNA (RIN={RIN})')
for ext in ['.pdf', '.png']:
        save(sample, f'./output/{DATE}_Bioanalyzer_{STRAIN}_total_RNA.{ext}')
# %%
