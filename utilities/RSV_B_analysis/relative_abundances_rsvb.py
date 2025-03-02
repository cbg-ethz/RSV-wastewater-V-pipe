import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker

df = pd.read_csv('../../preprint/data/lollipop/rsvb/deconvolved.csv', sep='\t')
viral_loads_2022_2023 = pd.read_csv("../../preprint/data/viral_load/normalized_viral_load/normalized_viral_loads.csv")


# Function to create and save a plot for each location
def plot_by_location(df, location, filename):
    # Filter data by location
    print(df)

    df_filtered = df[df['location'] == location]

    df_filtered = df_filtered.sort_values(by=['date'])
    df_filtered['date'] = pd.to_datetime(df_filtered['date'])

    # Set the plot style
    sns.set_theme(style="whitegrid")

    #palette = sns.color_palette("tab20", n_colors=df_filtered['variant'].nunique())
    variant_colors = {'B.D': 'olive',
                      'B.D.1': 'y',
                      'B.D.1.1': 'olivedrab',
                      'B.D.4': 'mediumorchid',
                      'B.D.4.1': 'plum',
                      'B.D.4.1.1': 'indigo',
                      'B.D.E.1': 'dodgerblue',
                      'B.D.E.2': 'powderblue',
                      'B.D.E.4': 'yellow',
                      'undetermined': 'black'}



    if location == "Zürich (ZH)":
        location = "Zurich (ZH)"
    elif location == "Genève (GE)":
        location = "Geneva (GE)"


    #variant_colors = dict(zip(df_filtered['variant'].unique(), palette))
    plt.figure(figsize=(12, 8))

    # Plot data for each variant
    for variant, group_data in df_filtered.groupby('variant'):
        #print(variant, group_data)
        plt.scatter(group_data['date'], group_data['proportion'], label=variant, color=variant_colors[variant])
        plt.plot(group_data['date'], group_data['proportion'], alpha=1, color=variant_colors[variant])
        plt.fill_between(group_data['date'], group_data['proportionLower'], group_data['proportionUpper'], alpha=0.1,
                         color=variant_colors[variant])


    # Customize the plot
    plt.xlabel('Date',fontsize=20)
    plt.ylabel('Proportion',fontsize=20)
    plt.title(f'2022-2023 - {location}', fontsize=25)
    plt.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=20,title_fontsize=20)
    #plt.tick_params(axis='both', which='major', labelsize=15)

    plt.gcf().subplots_adjust(bottom=0.2)  # Moves plot up to make space for x-axis label

    time= pd.to_datetime(np.unique(df_filtered[['date']]))
    plt.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=20,title_fontsize=20)
    plt.xticks(rotation=45)

    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=20)

    ax.set_xticks(time)
    ax.set_xticklabels(
        (time).strftime("%Y.%b.%d"),
        rotation=45,
        fontsize=20,
        ha='right')

    plt.tight_layout()

    # Save the plot
    plt.savefig(filename)
    # plt.show()


# Plot for "Genève (GE)"
plot_by_location(df, "Genève (GE)", "../../preprint/plots/relative_abundances/RSV_B/geneve_RSV_B_variants_over_time.pdf")

# Plot for "Zürich (ZH)"
plot_by_location(df, "Zürich (ZH)", "../../preprint/plots/relative_abundances/RSV_B/zurich_RSV_B_variants_over_time.pdf")


def plot_stacked_are_chart(df, location, treatment_plant, filename, viral_loads):

    data_rsv_viral_loads = viral_loads.loc[(viral_loads["wwtp"] == treatment_plant)].copy()
    data_rsv_viral_loads['Date'] = pd.to_datetime(data_rsv_viral_loads['sample_date'], format='%Y-%m-%d')
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


    df_filtered_new_stacked_area = df_filtered_new.mul(data_rsv_viral_loads['load_7days_median'], axis=0)

# Define custom colors for the plot
    custom_colors = [
            "olive",  # B.D
            "y",  # B.D.1
            "olivedrab",  # B.D.1.1
            "mediumorchid",  # B.D.4
            "plum",  # B.D.4.1
            "indigo",  # B.D.4.1.1
            "dodgerblue",  # B.D.E.1
            "powderblue",  # B.D.E.2
            "yellow", # B.D.E.4
            "black",  # Undetermined
    ]
    sns.set_theme(style="whitegrid")


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
# change to English names
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
        rotation=55,
        fontsize=20,
        ha='right')

    plt.xlabel('Date',fontsize=20)
    plt.ylabel('Flow-normalized viral load \n (gc/person/day)',fontsize=20)
    plt.title(f'2022-2023 {location}', fontsize=25)
    plt.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=25,title_fontsize=25)
    #plt.ylim(0, 3.5e07)

    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # Force scientific notation
    ax.tick_params(axis='y', labelsize=20)

    plt.tight_layout()

    plt.savefig(filename)


plot_stacked_are_chart(df=df,
                       location="Genève (GE)",
                       treatment_plant="STEP Aire",
                       filename="../../preprint/plots/relative_abundances/RSV_B/geneve_RSV_B_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2022_2023)

plot_stacked_are_chart(df=df,
                       location="Zürich (ZH)",
                       treatment_plant="ARA Werdhoelzli",
                       filename="../../preprint/plots/relative_abundances/RSV_B/zurich_RSV_B_variants_over_time_stacked_area.pdf",
                       viral_loads=viral_loads_2022_2023)

