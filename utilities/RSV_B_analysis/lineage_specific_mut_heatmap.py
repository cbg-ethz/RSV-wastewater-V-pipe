import re
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import json
import matplotlib.patches as mpatches


with open('../shared/RSV_data_analysis/rsv_definitions/RSVB_nucleotide_mutations_0.9.json', 'r') as file:
    clades_definitions = json.load(file)
clades = clades_definitions.keys()

swiss_clades = "B.D.E.1" # lineage for which all signature mutation frequencies are plotted

COVERAGE_THRESHOLD = 30

coverage = pd.read_csv("../../preprint/data/coverage/collected_rsv_coverage_rsv_b_2022_2023_PREPRINT.tsv")
timeline_tsv_mutation = pd.read_csv('../../preprint/data/timeline_tsv/RSV_B/timeline_mutation_rsv_b_2022_2023_PREPRINT.tsv',
                                    sep='\t', usecols=['submissionId',
                                                       'date',
                                                       'location',
                                                       'nucleotideMutationFrequency',
                                                       'aminoAcidMutationFrequency'])

timeline_tsv_mutation = timeline_tsv_mutation.sort_values(by=['location', 'date'],
                                                          ascending=[False, True])


timeline_tsv_mutation_sorted = {}
timeline_tsv_aminoacid_mutation_sorted = {}
for _, row in timeline_tsv_mutation.iterrows():#['date', 'nucleotideMutationFrequency']:
    nucleotide_mut_freq = row['nucleotideMutationFrequency']
    aminoacid_mut_freq = row['aminoAcidMutationFrequency']

    if pd.isna(nucleotide_mut_freq) or nucleotide_mut_freq == '':
      #  print(row)
        continue
    date = row['date']
    location = row['location']
    id = row['submissionId']

    sample_muts = json.loads(row['nucleotideMutationFrequency'])
    sample_aminoacid_muts = json.loads(row['aminoAcidMutationFrequency'])

    sample_muts_sorted = dict(sorted(sample_muts.items(),
                                               key=lambda x: int(re.findall(r'\d+', x[0])[0]) # sort by SNP location
                                               )
                                        )
    sample_aminoacid_muts_sorted = dict(sorted(sample_aminoacid_muts.items(),
                                               key=lambda x: int(re.findall(r'\d+', x[0])[-1]) #sort by the last number (SNP location)
                                               )
                                        )
    timeline_tsv_mutation_sorted[id] = sample_muts_sorted
    timeline_tsv_aminoacid_mutation_sorted[id] = sample_aminoacid_muts_sorted



timeline_tsv_mutation_sorted_df = pd.DataFrame.from_dict(timeline_tsv_mutation_sorted)
timeline_tsv_aminoacid_mutation_sorted_df = pd.DataFrame.from_dict(timeline_tsv_aminoacid_mutation_sorted)
df = timeline_tsv_mutation_sorted_df.transpose()

#

# drop deletions and insertions:
for mutation in df.columns:

    if (len(re.findall(r'[A-Za-z]+', mutation)[0]) > 1) or (len(re.findall(r'[A-Za-z]+', mutation)[1]) > 1):
        df = df.drop(columns=mutation)

# mutations that are not in B.D.E.1 definition are dropped -- we are plotting only lineage signature mutation frequencies
for mutation in df.columns:
    if mutation not in clades_definitions[swiss_clades]:
        #print(mutation)
        df = df.drop(columns=mutation)

#
print(df)
for mutation in clades_definitions[swiss_clades]:

    if mutation not in df.columns:
        df[mutation] = np.nan
        position = int(re.findall(r'\d+', mutation)[0])
        print(position)
        for sample in df.index:

            if (coverage.loc[(coverage['pos'] == position) & (coverage['sample'] == sample), 'coverage'].values[0] >= COVERAGE_THRESHOLD):
                if (pd.isna(df.loc[sample, mutation])):
                    df.loc[sample, mutation] = 0.0

        # if coverage at the position is below the COVERAGE_THRESHOLD -> missing value
            else:
                df.loc[sample, mutation] = np.nan


df.index = [re.sub('_3x', '', str(ind)) for ind in df.index]
df.index = [re.sub(r'^16_', 'GE_', str(ind)) for ind in df.index]
df.index = [re.sub(r'^10_', 'ZH_', str(ind)) for ind in df.index]

# drop mutations in non-coding regions and synonymous mutations (to save the place in full-genome plot):
df_full_genome_plot = df

# dataframe of mutation definitions
sorted_columns = sorted(df_full_genome_plot.columns, key = lambda x: int(re.findall(r'\d+', x)[0]))
df_full_genome_plot = df_full_genome_plot[sorted_columns]

mutations_df = pd.DataFrame(0, index=[swiss_clades], columns=df_full_genome_plot.columns)
for clade in mutations_df.index:
    for mut in mutations_df.columns:
        mutations_df.at[clade, mut] = 1 if mut in clades_definitions[clade] else 0


sns.set(rc={'figure.figsize': (30, 15)})
sns.set_style("white")
plt.grid(True, linewidth=0.1, color='gray')

fig, axs = plt.subplots(figsize=(40, 20))

sns.heatmap(df_full_genome_plot, ax=axs, yticklabels=df_full_genome_plot.index.to_list(),
            linecolor="black",
            linewidths=0.0,
            cmap=sns.color_palette("Blues", as_cmap=True),
            cbar_kws={"shrink": 0.5, "aspect": 10,"label": "Frequency", "pad":0.01})


# Access the colorbar object
colorbar = axs.collections[0].colorbar

# Customize the tick labels and label font size
colorbar.ax.tick_params(labelsize=28)  # Adjust tick label font size
colorbar.set_label("Frequency", fontsize=32)  # Adjust the colorbar label size

# Iterate over the X-axis ticks and color based on "GE" or "ZH"
for ytick in axs.get_yticklabels():
    label_text = ytick.get_text()
    if "GE" in label_text:
        #ytick.set_color("green")
        ytick.set_backgroundcolor("lightgreen")  # Background color

    elif "ZH" in label_text:
        #ytick.set_color("blue")
        ytick.set_backgroundcolor("lightblue")  # Background color


axs.set_facecolor("#ffe6e6")
axs.set_xticks([x + 0.5 for x in range(df_full_genome_plot.shape[1])])
axs.set_xticklabels([idx.split('_')[0] for idx in df_full_genome_plot.transpose().index],
                    fontsize=25,
                    rotation=90,
                    ha='center',
                    va='top')
axs.tick_params(axis='y', labelsize=35)

plt.tight_layout()

color_ranges = [
    (57, 476, "red"),
    (584, 958, "khaki"),
    (1097, 2272, "green"),
    (2305, 3030, "purple"),
    (3154, 3990, "coral"),
    (4259, 4456, "cyan"),
    (4646, 5578, "peru"),
    (5676, 7400, "darkcyan"),
    (7627, 8452, "orange"),
    (8532, 15029, "blue"),
]
for xtick in axs.get_xticklabels():
    #xtick.set_fontweight("bold")
    xtick_text = xtick.get_text().split('_')[0]
    uppercase_letters = (re.findall(r'[A-Z]', xtick_text))  # Extract genome position

    if len(uppercase_letters) >= 2:  # Ensure there are at least two letters to compare
       # print(uppercase_letters)
        aa_ref = uppercase_letters[0]  # Reference genome position
        aa_alt = uppercase_letters[1]  # Alternative genome position

        if aa_ref != aa_alt:
            xtick.set_color('black')
            xtick.set_fontweight("bold")



for xtick in axs.get_xticklabels():
    pos_in_genome = int(re.findall(r'\d+', xtick.get_text())[0])  # Extract genome position
    for start, end, color in color_ranges:
        if start < pos_in_genome < end:
            xtick.set_color(color)

            xtick.set_fontweight("bold")
            break

# Define the color-coded regions as a list of tuples
color_regions = [
    ("NS1 (57-476)", "red"),
    ("NS2 (584-958)", "khaki"),
    ("N (1097-2272)", "green"),
    ("P (2305-3030)", "purple"),
    ("M (3154-3990)", "coral"),
    ("SH (4259-4456)", "cyan"),
    ("G (4646-5578)", "peru"),
    ("F (5676-7400)", "darkcyan"),
    ("M2-1, M2-2 (7627-8452)", "orange"),
    ("L (8518-15018)", "blue"),
]

# Create a list of mpatches.Patch objects for the legend
legend_patches = [mpatches.Patch(color=color, label=label) for label, color in color_regions]

# Add the legend to the plot (you can customize location and other properties)
axs.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.08),
              ncol=5, fontsize=35, frameon=False, title="Genome Regions",title_fontsize=35)

plt.suptitle("RSV-B B.D.E.1 lineage, 2022-2023 season",
             fontsize=30, fontweight='bold', y=1.02)

fig.savefig("../../preprint/plots/mutation_heatmap/BDE1_22_23_final.pdf", format="pdf", bbox_inches="tight")