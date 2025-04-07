import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data
batch = '20250307_2418539124_and_20250307_2418653583_experiment'
subbatch ='20250307_2418653583'
reference = "EPI_ISL_1653999"
subtype='B'

#reference = "EPI_ISL_1653999"
#subtype='B'

#df = pd.read_csv(f"../../data/all_data/{batch}/collected_rsv_coverage_regular_monitoring_rsv_{subtype}_2024_2025_{batch}_{reference}.tsv")
df =pd.read_csv(f"/Users/arimaite/Documents/GitHub/RSV_wastewater/regular_monitoring_2024_2025/data/all_data/{batch}/collected_rsv_coverage_{subbatch}_{reference}.tsv")
#print(df)
# Calculate log coverage and replace -inf with 0
df['log_coverage'] = np.log10(df['coverage']).replace(-np.inf, 0)

df = df.drop(columns='coverage')
print(df)
# Filter for samples that start with '10', '16'
lugano_05 = df[df['sample'].str.contains('_05_2025_')]
zurich_10 = df[df['sample'].str.contains('_10_2025_')]
basel_15 = df[df['sample'].str.contains('_15_2025_')]
geneva_16 = df[df['sample'].str.contains('_16_2025_')]
chur_17 = df[df['sample'].str.contains('_17_2025_')]
laupen_25 = df[df['sample'].str.contains('_25_2025_')]
# Group the data by sample
grouped_lugano_05 = lugano_05.groupby('sample')
grouped_zurich_10 = zurich_10.groupby('sample')
grouped_basel_15 = basel_15.groupby('sample')
grouped_geneva_16 = geneva_16.groupby('sample')
grouped_chur_17 = chur_17.groupby('sample')
grouped_laupen_25 = laupen_25.groupby('sample')
print(geneva_16)
for sample_name, group in grouped_geneva_16:
    print(f"Sample: {sample_name}")
    print(group.head())  # Display the first few rows of each group

############


# ZURICH
# Number of samples to plot
num_samples = len(grouped_zurich_10)

# Set number of columns and rows dynamically
num_cols = 1  # Fixed number of columns
num_rows = 6  # Calculate rows based on number of samples

# Create a figure with subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(20,15))
axes = axes.flatten()

# Create a plot for each sample in a subplot
for ax, (sample_name, group) in zip(axes, grouped_zurich_10):
    ax.plot(group['pos'], group['log_coverage'], linestyle='-', linewidth=0.5)
    ax.fill_between(group['pos'], group['log_coverage'], color='lightblue', alpha=0.5)

    ax.set_title(f'{sample_name}', fontsize=25, fontweight='bold')  # Larger title, bold
    ax.set_xlabel('Position', fontsize=18, fontweight='bold')  # Bold and larger x-axis label
    ax.set_ylabel('Log Coverage', fontsize=18, fontweight='bold')  # Bold and larger y-axis label
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=16, width=2, length=10)  # Major ticks
    #ax.set_ylim(-1, 6)  # Set appropriate limits based on your data

    # Set thicker axes spines
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

# Turn off any unused subplots
for ax in axes[num_samples:]:
    ax.axis('off')  # Hide unused axes

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the figure with all subplots
plt.savefig(f'../../plots/all_data/{batch}/rsv_{subtype}_log_coverage_plots_zurich_{reference}_{subbatch}.png')  # Save as PNG or other formats
plt.close()  # Close the figure to free memory


# GENEVA
# Number of samples to plot
num_samples = len(grouped_geneva_16)

# Set number of columns and rows dynamically
num_cols = 1  # Fixed number of columns
num_rows = 6  # Calculate rows based on number of samples

# Create a figure with subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(20,15))
axes = axes.flatten()

# Create a plot for each sample in a subplot
for ax, (sample_name, group) in zip(axes, grouped_geneva_16):
    ax.plot(group['pos'], group['log_coverage'], linestyle='-', linewidth=0.5)
    ax.fill_between(group['pos'], group['log_coverage'], color='lightblue', alpha=0.5)

    ax.set_title(f'{sample_name}', fontsize=25, fontweight='bold')  # Larger title, bold
    ax.set_xlabel('Position', fontsize=18, fontweight='bold')  # Bold and larger x-axis label
    ax.set_ylabel('Log Coverage', fontsize=18, fontweight='bold')  # Bold and larger y-axis label
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=16, width=2, length=10)  # Major ticks
    # Set y-axis to logarithmic scale

    # Set the logarithmic axis formatter for y-axis to show 10^x format

    # Set thicker axes spines
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

# Turn off any unused subplots
for ax in axes[num_samples:]:
    ax.axis('off')  # Hide unused axes

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the figure with all subplots
plt.savefig(f'../../plots/all_data/{batch}/rsv_{subtype}_log_coverage_plots_geneva_{reference}_{subbatch}.png')  # Save as PNG or other formats
plt.close()  # Close the figure to free memory




# LUGANO
# Number of samples to plot
num_samples = len(grouped_lugano_05)

# Set number of columns and rows dynamically
num_cols = 1  # Fixed number of columns
num_rows = 6  # Calculate rows based on number of samples

# Create a figure with subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(20,15))
axes = axes.flatten()

# Create a plot for each sample in a subplot
for ax, (sample_name, group) in zip(axes, grouped_lugano_05):
    ax.plot(group['pos'], group['log_coverage'], linestyle='-', linewidth=0.5)
    ax.fill_between(group['pos'], group['log_coverage'], color='lightblue', alpha=0.5)

    ax.set_title(f'{sample_name}', fontsize=25, fontweight='bold')  # Larger title, bold
    ax.set_xlabel('Position', fontsize=18, fontweight='bold')  # Bold and larger x-axis label
    ax.set_ylabel('Log Coverage', fontsize=18, fontweight='bold')  # Bold and larger y-axis label
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=16, width=2, length=10)  # Major ticks
    #ax.set_ylim(-1, 6)  # Set appropriate limits based on your data

    # Set thicker axes spines
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

# Turn off any unused subplots
for ax in axes[num_samples:]:
    ax.axis('off')  # Hide unused axes

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the figure with all subplots
plt.savefig(f'../../plots/all_data/{batch}/rsv_{subtype}_log_coverage_plots_lugano_{reference}_{subbatch}.png')  # Save as PNG or other formats
plt.close()  # Close the figure to free memory




# BASEL
# Number of samples to plot
num_samples = len(grouped_basel_15)

# Set number of columns and rows dynamically
num_cols = 1  # Fixed number of columns
num_rows = 6  # Calculate rows based on number of samples

# Create a figure with subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(20,15))
axes = axes.flatten()

# Create a plot for each sample in a subplot
for ax, (sample_name, group) in zip(axes, grouped_basel_15):
    ax.plot(group['pos'], group['log_coverage'], linestyle='-', linewidth=0.5)
    ax.fill_between(group['pos'], group['log_coverage'], color='lightblue', alpha=0.5)

    ax.set_title(f'{sample_name}', fontsize=25, fontweight='bold')  # Larger title, bold
    ax.set_xlabel('Position', fontsize=18, fontweight='bold')  # Bold and larger x-axis label
    ax.set_ylabel('Log Coverage', fontsize=18, fontweight='bold')  # Bold and larger y-axis label
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=16, width=2, length=10)  # Major ticks
    #ax.set_ylim(-1, 6)  # Set appropriate limits based on your data

    # Set thicker axes spines
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

# Turn off any unused subplots
for ax in axes[num_samples:]:
    ax.axis('off')  # Hide unused axes

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the figure with all subplots
plt.savefig(f'../../plots/all_data/{batch}/rsv_{subtype}_log_coverage_plots_basel_{reference}_{subbatch}.png')  # Save as PNG or other formats
plt.close()  # Close the figure to free memory





# CHUR
# Number of samples to plot
num_samples = len(grouped_chur_17)

# Set number of columns and rows dynamically
num_cols = 1  # Fixed number of columns
num_rows = 6  # Calculate rows based on number of samples

# Create a figure with subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(20,15))
axes = axes.flatten()

# Create a plot for each sample in a subplot
for ax, (sample_name, group) in zip(axes, grouped_chur_17):
    ax.plot(group['pos'], group['log_coverage'], linestyle='-', linewidth=0.5)
    ax.fill_between(group['pos'], group['log_coverage'], color='lightblue', alpha=0.5)

    ax.set_title(f'{sample_name}', fontsize=25, fontweight='bold')  # Larger title, bold
    ax.set_xlabel('Position', fontsize=18, fontweight='bold')  # Bold and larger x-axis label
    ax.set_ylabel('Log Coverage', fontsize=18, fontweight='bold')  # Bold and larger y-axis label
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=16, width=2, length=10)  # Major ticks
    #ax.set_ylim(-1, 6)  # Set appropriate limits based on your data

    # Set thicker axes spines
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

# Turn off any unused subplots
for ax in axes[num_samples:]:
    ax.axis('off')  # Hide unused axes

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the figure with all subplots
plt.savefig(f'../../plots/all_data/{batch}/rsv_{subtype}_log_coverage_plots_chur_{reference}_{subbatch}.png')  # Save as PNG or other formats
plt.close()  # Close the figure to free memory


# LAUPEN
# Number of samples to plot
num_samples = len(grouped_laupen_25)

# Set number of columns and rows dynamically
num_cols = 1  # Fixed number of columns
num_rows = 6  # Calculate rows based on number of samples

# Create a figure with subplots
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(20,15))
axes = axes.flatten()

# Create a plot for each sample in a subplot
for ax, (sample_name, group) in zip(axes, grouped_laupen_25):
    ax.plot(group['pos'], group['log_coverage'], linestyle='-', linewidth=0.5)
    ax.fill_between(group['pos'], group['log_coverage'], color='lightblue', alpha=0.5)

    ax.set_title(f'{sample_name}', fontsize=25, fontweight='bold')  # Larger title, bold
    ax.set_xlabel('Position', fontsize=18, fontweight='bold')  # Bold and larger x-axis label
    ax.set_ylabel('Log Coverage', fontsize=18, fontweight='bold')  # Bold and larger y-axis label
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=16, width=2, length=10)  # Major ticks
    #ax.set_ylim(-1, 6)  # Set appropriate limits based on your data

    # Set thicker axes spines
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

# Turn off any unused subplots
for ax in axes[num_samples:]:
    ax.axis('off')  # Hide unused axes

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the figure with all subplots
plt.savefig(f'../../plots/all_data/{batch}/rsv_{subtype}_log_coverage_plots_laupen_{reference}_{subbatch}.png')  # Save as PNG or other formats
plt.close()  # Close the figure to free memory