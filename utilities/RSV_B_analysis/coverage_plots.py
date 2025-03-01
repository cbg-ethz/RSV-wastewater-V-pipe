import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# Load the data
df = pd.read_csv("../../preprint/data/coverage/collected_rsv_coverage_rsv_b_2022_2023_PREPRINT.tsv")

print(np.unique(df['sample']))
print(df)
# add one to avoid calculating log of zero
df['log_coverage'] = np.log10(df['coverage']+1)

df = df.drop(columns='coverage')
print(df)
# Filter for samples that start with '10', '16' (code for the corresponding locations)
df_zurich = df[df['sample'].str.startswith('10')]
df_geneva = df[df['sample'].str.startswith('16')]

# Group the data by sample
grouped_zh = df_zurich.groupby('sample')
grouped_ge = df_geneva.groupby('sample')


# ZURICH
# Number of samples to plot
num_samples = len(grouped_zh)

num_cols = 4  # number of columns
num_rows = 5  # number of rows

# Create a figure with subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(35, 20))
axes = axes.flatten()

# Create a plot for each sample in a subplot
for ax, (sample_name, group) in zip(axes, grouped_zh):

    sample_name = re.sub(r'^10', 'ZH', sample_name)
    #print(sample_name)
    sample_name = re.sub('_3x', '', sample_name)
    #print(sample_name)


    ax.plot(group['pos'], group['log_coverage'], linestyle='-', linewidth=0.5)
    ax.fill_between(group['pos'], group['log_coverage'], color='lightblue', alpha=0.5)

    ax.set_title(f'{sample_name}', fontsize=25, fontweight='bold')
    ax.set_xlabel('Position', fontsize=25, fontweight='bold')
    ax.set_ylabel('Log Coverage', fontsize=25, fontweight='bold')
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=20, width=2, length=10)
    ax.set_ylim(-1, 7)

    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

# Turn off any unused subplots
for ax in axes[num_samples:]:
    ax.axis('off')

# set tight layout
plt.tight_layout()

# Save the figure with all subplots
plt.savefig('../../preprint/plots/coverage/log_coverage_plots_zurich_2022_2023.png')  # Save as PNG or other formats
plt.close()  # Close the figure


# GENEVA
# Number of samples to plot
num_samples = len(grouped_ge)

num_cols = 4  # number of columns
num_rows = 5  # number of rows

# Create a figure with subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(35, 20))
axes = axes.flatten()

# Create a plot for each sample in a subplot
for ax, (sample_name, group) in zip(axes, grouped_ge):
    sample_name = re.sub(r'^16', 'GE', sample_name)


    ax.plot(group['pos'], group['log_coverage'], linestyle='-', linewidth=0.5)
    ax.fill_between(group['pos'], group['log_coverage'], color='lightblue', alpha=0.5)

    ax.set_title(f'{sample_name}', fontsize=25, fontweight='bold')
    ax.set_xlabel('Position', fontsize=25, fontweight='bold')
    ax.set_ylabel('Log Coverage', fontsize=25, fontweight='bold')
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=20, width=2, length=10)

    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

# Turn off any unused subplots
for ax in axes[num_samples:]:
    ax.axis('off')

# set tight layout
plt.tight_layout()

# Save the figure with all subplots
plt.savefig('../../preprint/plots/coverage/log_coverage_plots_geneva_2022_2023.png')  # Save as PNG or other formats
plt.close()  # Close the figure
