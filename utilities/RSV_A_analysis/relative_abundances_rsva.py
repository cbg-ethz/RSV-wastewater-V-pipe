import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker

df = pd.read_csv('../../preprint/data/lollipop/rsva/deconvolved.csv', sep='\t')
viral_loads_2023_2024 = pd.read_csv("../../preprint/data/viral_load/normalized_viral_load/2023-11-01-2024-03-01WISE-dPCR.csv")

# Function to create and save a plot for each location
def plot_by_location(df, location, filename):
    # Filter data by location

    df_filtered = df[df['location'] == location]

    df_filtered = df_filtered.sort_values(by=['date'])
    df_filtered['date'] = pd.to_datetime(df_filtered['date'])

    # Set the plot style
    sns.set_theme(style="whitegrid")
    #palette = sns.color_palette("tab20", n_colors=df_filtered['variant'].nunique())

    variant_colors = {'A.D.1': 'olive',
                      'A.D.1.5': 'y',
                      'A.D.1.6': 'olivedrab',
                      'A.D.2.1': 'darkcyan',
                      'A.D.3': 'mediumorchid',
                      'A.D.3.1': 'plum',
                      'A.D.5.1': 'dodgerblue',
                      'A.D.5.2': 'powderblue',
                      'undetermined': 'black'}

    print(variant_colors)
    plt.figure(figsize=(12, 8))

    if location == "Zürich (ZH)":
        location = "Zurich (ZH)"
    elif location == "Genève (GE)":
        location = "Geneva (GE)"


    # Plot data for each variant
    for variant, group_data in df_filtered.groupby('variant'):
        #print(variant, group_data)
        plt.scatter(group_data['date'], group_data['proportion'], label=variant,color=variant_colors[variant])
        plt.plot(group_data['date'], group_data['proportion'], alpha=1, color=variant_colors[variant])
        plt.fill_between(group_data['date'], group_data['proportionLower'], group_data['proportionUpper'], alpha=0.1, color=variant_colors[variant])

    # Customize the plot
    plt.xlabel('Date',fontsize=20)
    plt.ylabel('Proportion',fontsize=20)
    plt.title(f'2023-2024 - {location}', fontsize=25)
    plt.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=20,title_fontsize=20)
    #plt.tick_params(axis='both', which='major', labelsize=15)

    plt.gcf().subplots_adjust(bottom=0.2)  # Moves plot up to make space for x-axis label

    plt.xticks(rotation=45)
    ax = plt.gca()  # Get the current Axes
    ax.set_xticks(df_filtered['date'])
    ax.set_xticklabels(
        df_filtered['date'].dt.strftime("%Y.%b.%d"),  # Ensure 'date' is a datetime type
        rotation=55,
        fontsize=20,
        ha='right')
    plt.tight_layout()

    ax.tick_params(axis='y', labelsize=20)

    # Save the plot
    plt.savefig(filename)


# Plot for "Genève (GE)"
plot_by_location(df, "Genève (GE)", "../../preprint/plots/relative_abundances/RSV_A/rsva_geneve_variants_over_time.pdf")

# Plot for "Zürich (ZH)"
plot_by_location(df, "Zürich (ZH)", "../../preprint/plots/relative_abundances/RSV_A/rsva_zurich_variants_over_time.pdf")

def plot_stacked_area_chart(df, location, treatment_plant, filename, viral_loads):

    data_rsv_viral_loads = viral_loads.loc[((viral_loads["Virus"] == "RSV-N")|(viral_loads["Virus"] == "RSV-M")) &
                                           (viral_loads["Treatment Plant"] == treatment_plant)].copy()
    data_rsv_viral_loads['Date'] = pd.to_datetime(data_rsv_viral_loads['Date'], format='%Y-%m-%d')
    data_rsv_viral_loads.sort_values(by='Date', inplace=True)


    df_filtered = df[df['location'] == location]
    time= pd.to_datetime(np.unique(df_filtered[['date']]))

    df_filtered_new = pd.DataFrame()

    for variant in np.unique(df_filtered['variant']):
        new_col = df_filtered.loc[df_filtered['variant'] == variant, ['proportion','date']]
        df_filtered_new[variant] = new_col.reset_index()['proportion']

    df_filtered_new.index = time

    data_rsv_viral_loads = data_rsv_viral_loads.loc[data_rsv_viral_loads['Date'].isin(time)].reset_index()
    data_rsv_viral_loads.index = df_filtered_new.index


    df_filtered_new_stacked_area = df_filtered_new.mul(data_rsv_viral_loads['7-day Median Viral Load (gc/person/day)'], axis=0)

# Define custom colors for the plot
    custom_colors = [
        "olive",  # A.D.1
        "y",  # A.D.1.5
        "olivedrab",  # A.D.1.6
        "darkcyan",  # A.D.2.1
        "mediumorchid",  # A.D.3
        "plum",  # A.D.3.1
        "dodgerblue",  # A.D.5.1
        "powderblue",  # A.D.5.2
        "black",  # Undetermined
    ]
    sns.set_theme(style="whitegrid")

    # Prepare data for stack plotting
   # areas = np.stack(df_filtered_new.values, axis=-1)

    # To get plots of relative abundances -> uncomment the next line:
    areas = np.stack(df_filtered_new_stacked_area.values, axis=-1)
# Create the figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))  # Adjust figure size for better readability
# Plot the stackplot
    ax.stackplot(
        pd.to_datetime(df_filtered_new_stacked_area.index),  # x-axis values
        areas,                   # y-axis areas
        labels=np.unique(df_filtered['variant']),  # Labels for each stack
        colors=custom_colors,    # Custom colors
        alpha=0.7                # Slight transparency for better visualization
    )

# change names to english
    if location == "Zürich (ZH)":
        location = "Zurich (ZH)"
    elif location == "Genève (GE)":
        location = "Geneva (GE)"


    #ax.tick_params(axis='x', rotation=45)  # Rotate x-axis labels for clarity
    #ax.tick_params(axis='both', which='major', labelsize=15)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=10, title="Variants")
    plt.xticks(rotation=45)

    ax = plt.gca()
    ax.set_xticks(time)
    ax.set_xticklabels(
        (time).strftime("%Y.%b.%d"),  # Ensure 'date' is a datetime type
        rotation=45,
        fontsize=20,
        ha='right')

    plt.xlabel('Date',fontsize=20)
    plt.ylabel('Flow-normalized viral load \n (gc/person/day)',fontsize=20)
    plt.title(f'2023-2024 {location}', fontsize=25)
    plt.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=20,title_fontsize=20)
    #plt.tick_params(axis='both', which='major', labelsize=15)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # Force scientific notation
    ax.tick_params(axis='y', labelsize=20)

    plt.tight_layout()

    plt.savefig(filename)

plot_stacked_area_chart(df=df,
                       location="Genève (GE)",
                       treatment_plant ="STEP Aire",
                       filename="../../preprint/plots/relative_abundances/RSV_A/rsva_geneve_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2023_2024)

plot_stacked_area_chart(df=df,
                       location="Zürich (ZH)",
                       treatment_plant ="ARA Werdhoelzli",
                       filename="../../preprint/plots/relative_abundances/RSV_A/rsva_zurich_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2023_2024)

