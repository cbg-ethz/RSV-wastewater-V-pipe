import re
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import argparse
import requests
import json
import matplotlib.patches as mpatches


with open('/Users/arimaite/Documents/GitHub/RSV_wastewater/utilities/shared/RSV_data_analysis/rsv_definitions/RSVA_nucleotide_mutations_0.9.json', 'r') as file:
    clades_definitions = json.load(file)
clades = clades_definitions.keys()

swiss_clades = ["A.D.1", "A.D.1.5", "A.D.1.6", "A.D.2.1", "A.D.3", "A.D.3.1", "A.D.5.1", "A.D.5.2"]

timeline_tsv_mutation = pd.read_csv('../../data/all_data/20250321_2429695737_regular_monitoring/20250321_2429695737_EPI_ISL_412866_Mutations_Dashboard.tsv',
                                    sep='\t',
                                    usecols=['date',
                                             'location',
                                             'nucleotideMutationFrequency',
                                             'aminoAcidMutationFrequency'])
timeline_tsv_mutation = timeline_tsv_mutation.sort_values(by=['location', 'date'],
                                                          ascending=[False, True])


timeline_tsv_mutation_sorted = {}
timeline_tsv_aminoacid_mutation_sorted = {}
for _, row in timeline_tsv_mutation.iterrows():
    nucleotide_mut_freq = row['nucleotideMutationFrequency']
    aminoacid_mut_freq = row['aminoAcidMutationFrequency']

    if pd.isna(nucleotide_mut_freq) or nucleotide_mut_freq == '':

        continue
    date = row['date']
    location = row['location']
# Mutations are already sorted by nt position by TSV generating script
    sample_muts_sorted = json.loads(row['nucleotideMutationFrequency'])
    sample_aminoacid_muts_sorted = json.loads(row['aminoAcidMutationFrequency'])


    timeline_tsv_mutation_sorted[str(date)+"_"+location] = sample_muts_sorted
    timeline_tsv_aminoacid_mutation_sorted[str(date)+"_"+location] = sample_aminoacid_muts_sorted


timeline_tsv_mutation_sorted_df = pd.DataFrame.from_dict(timeline_tsv_mutation_sorted)
timeline_tsv_aminoacid_mutation_sorted_df = pd.DataFrame.from_dict(timeline_tsv_aminoacid_mutation_sorted)
df = timeline_tsv_mutation_sorted_df.transpose()
df_aa = timeline_tsv_aminoacid_mutation_sorted_df.transpose()

# # drop deletions and insertions:
# for mutation in df.columns:
#
#     if (len(re.findall(r'[A-Za-z]+', mutation)[0]) > 1) or (len(re.findall(r'[A-Za-z]+', mutation)[1]) > 1):
#         df = df.drop(columns=mutation)
#         df_aa = df_aa.loc[:, ~df_aa.columns.str.contains(f'_{re.escape(mutation)}')]


# dataframe of lineages definitions
mutations_df = pd.DataFrame(0, index=swiss_clades, columns=df.columns)
for clade in mutations_df.index:
    for mut in mutations_df.columns:
        mutations_df.at[clade, mut] = 1 if mut in clades_definitions[clade] else 0


sns.set_style("white")

plt.grid(True, linewidth=0.1, color='gray')

sns.set(rc={'figure.figsize': (40, 20)})

sns.set_style("white")

plt.grid(True, linewidth=0.1, color='gray')

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(40, 20), gridspec_kw={'height_ratios': [8, 2]})
df = df.apply(pd.to_numeric)

sns.heatmap(df, ax=axs[0], yticklabels=df.index.to_list(),linecolor="black", linewidths=0.0,
            cmap=sns.color_palette("Blues", as_cmap=True),cbar_kws={"shrink": 0.5, "aspect": 10,"label": "Frequency","pad":0.01})


# Access the colorbar object
colorbar = axs[0].collections[0].colorbar

# Customize the tick labels and label font size
colorbar.ax.tick_params(labelsize=20)  # Adjust tick label font size
colorbar.set_label("Frequency", fontsize=20)  # Adjust the colorbar label size



sns.heatmap(mutations_df, ax=axs[1], yticklabels=swiss_clades, linecolor="black", linewidths=0.0,
            cmap=sns.color_palette("Blues", as_cmap=True), cbar_kws={"shrink": 0.5, "pad": 0.01,"aspect": 10})
# make colorbar not visible in the second plot
colorbar2 = axs[1].collections[0].colorbar
colorbar2.ax.set_visible(False)

# Iterate over the X-axis ticks and color based on location

for ytick in axs[0].get_yticklabels():
    label_text = ytick.get_text()
    #ytick.set_fontweight("bold")
    if "Zurich" in label_text:
        ytick.set_color("green")
    elif "Lugano" in label_text:
        ytick.set_color("purple")
    elif "Laupen" in label_text:
        ytick.set_color("cyan")
    elif "Geneva" in label_text:
        ytick.set_color("blue")
    elif "Chur" in label_text:
        ytick.set_color("red")
    elif "Basel" in label_text:
        ytick.set_color("orange")


axs[0].set_facecolor("#ffe6e6")
axs[0].set_xticks([x + 0.5 for x in range(df.shape[1])])
axs[0].set_xticklabels([idx.split('_')[0] for idx in df_aa.transpose().index], fontsize=15, rotation=90, ha='center', va='top')
axs[0].tick_params(axis='y', labelsize=10)

plt.tight_layout()

bold_regions = [
    (4652, 5617, "G"),  # G region
    (5697, 7421, "F"),  # F region
]

axs[1].set_xticks([x + 0.5 for x in range(df.shape[1])])
axs[1].set_xticklabels(df.transpose().index, fontsize=15, rotation=90, ha='center', va='top')
axs[1].set_facecolor("#ffe6e6")
axs[1].tick_params(axis='y', labelsize=32)
color_ranges = [
    (70, 489, "red"),
    (599, 973, "khaki"),
    (1111, 2286, "green"),
    (2318, 3043, "purple"),
    (3226, 3996, "coral"),
    (4266, 4460, "cyan"),
    (4652, 5617, "peru"),
    (5697, 7421, "darkcyan"),
    (7640, 8465, "orange"),
    (8532, 15029, "blue"),
]


# Iterate over the X-axis ticks and assign colors based on ranges
for index, xtick in enumerate(axs[1].get_xticklabels()):
    pos_in_genome = int(re.findall(r'\d+', xtick.get_text())[0])  # Extract genome position
    for start, end, color in color_ranges:
        if start < pos_in_genome < end:
            xtick.set_color(color)
            xtick.set_fontweight("bold")
            axs[0].get_xticklabels()[index].set_color(color)
            break

# Define the color-coded regions as a list of tuples
color_regions = [
    ("NS1 (70–489)", "red"),
    ("NS2 (599–973)", "khaki"),
    ("N (1111–2286)", "green"),
    ("P (2318–3043)", "purple"),
    ("M (3226–3996)", "coral"),
    ("SH (4266–4460)", "cyan"),
    ("G (4652–5617)", "peru"),
    ("F (5697–7421)", "darkcyan"),
    ("M2-1, M2-2 (7640-8465)", "orange"),
    ("L (8532–15029)", "blue"),
]

# Create a list of mpatches.Patch objects for the legend
legend_patches = [mpatches.Patch(color=color, label=label) for label, color in color_regions]

# Add the legend to the plot (you can customize location and other properties)
axs[1].legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.2),
              ncol=4, fontsize=30, frameon=False, title="Genome Regions", title_fontsize=30)
#plt.subplots_adjust(hspace=0.0)  # Adjust the space as needed
plt.suptitle("Mutation frequencies (RSV-A, 2024-2025 season)",
             fontsize=30, fontweight='bold', y=1.02)

fig.savefig("../../plots/all_data/20250321_2429695737_regular_monitoring/nexus_20250321_2429695737_full_genome_rsva_nonsyn.pdf", format="pdf", bbox_inches="tight")

print(df_aa)
df_F_gene = df.loc[:, [col for col in df.columns if 5697 < int(re.findall(r'\d+', col)[0]) < 7421]]
df_aa_F_gene = df_aa.loc[:, [col for col in df_aa.columns if col.split(":")[0]=="F"]]
print(df_aa_F_gene)


mutations_df_F_gene = pd.DataFrame(0, index=swiss_clades, columns=df_F_gene.columns)
for clade in mutations_df_F_gene.index:
    for mut in mutations_df_F_gene.columns:
        mutations_df_F_gene.at[clade, mut] = 1 if mut in clades_definitions[clade] else 0

sns.set_style("white")
plt.grid(True, linewidth=0.1, color='gray')
#sns.set(rc={'figure.figsize': (30, 20)})
sns.set_style("white")
plt.grid(True, linewidth=0.1, color='gray')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(20, 30), gridspec_kw={'height_ratios': [7, 1]})

sns.heatmap(df_F_gene, ax=axs[0], yticklabels=df_F_gene.index.to_list(),linecolor="black", linewidths=0.0,
            cmap=sns.color_palette("Blues", as_cmap=True),cbar_kws={"shrink": 0.5, "aspect": 10,"pad": 0.02,"label": "Frequency"})

# Access the colorbar object
colorbar = axs[0].collections[0].colorbar

# Customize the tick labels and label font size
colorbar.ax.tick_params(labelsize=25)  # Adjust tick label font size
colorbar.set_label("Frequency", fontsize=35)  # Adjust the colorbar label size


sns.heatmap(mutations_df_F_gene, ax=axs[1], yticklabels=swiss_clades, linecolor="black", linewidths=0.0,
            cmap=sns.color_palette("Blues", as_cmap=True), cbar_kws={"shrink": 0.5, "aspect": 10,"pad": 0.02})
# make colorbar not visible in the second plot
colorbar2 = axs[1].collections[0].colorbar
colorbar2.ax.set_visible(False)
# Iterate over the X-axis ticks and color based on location
for ytick in axs[0].get_yticklabels():
    label_text = ytick.get_text()
    #ytick.set_fontweight("bold")
    if "Zurich" in label_text:
        ytick.set_color("green")
    elif "Lugano" in label_text:
        ytick.set_color("purple")
    elif "Laupen" in label_text:
        ytick.set_color("cyan")
    elif "Geneva" in label_text:
        ytick.set_color("blue")
    elif "Chur" in label_text:
        ytick.set_color("red")
    elif "Basel" in label_text:
        ytick.set_color("orange")

axs[0].set_facecolor("#ffe6e6")
axs[0].set_xticks([x + 0.5 for x in range(df_F_gene.shape[1])])
axs[0].set_xticklabels(df_aa_F_gene.transpose().index, fontsize=30, rotation=90, ha='center', va='top')
axs[0].tick_params(axis='y', labelsize=15)

plt.tight_layout()


#plt.show()
bold_regions = [
    (5697, 7421, "F")  # F region
]

axs[1].set_xticks([x + 0.5 for x in range(df_F_gene.shape[1])])
axs[1].set_xticklabels(df_F_gene.transpose().index, fontsize=30, rotation=90, ha='center', va='top')
axs[1].set_facecolor("#ffe6e6")
axs[1].tick_params(axis='y', labelsize=30, rotation=0)

color_ranges = [
    (5697, 7421, "darkcyan")
]

# Iterate over the X-axis ticks and assign colors based on ranges
for index, xtick in enumerate(axs[1].get_xticklabels()):
    pos_in_genome = int(re.findall(r'\d+', xtick.get_text())[0])  # Extract genome position
    for start, end, color in color_ranges:
        if start < pos_in_genome < end:
            xtick.set_color(color)
            xtick.set_fontweight("bold")
            axs[0].get_xticklabels()[index].set_fontweight("bold")
            break


# Define the color-coded regions as a list of tuples
color_regions = [
    ("F (5697–7421)", "darkcyan")
]

# Create a list of mpatches.Patch objects for the legend
legend_patches = [mpatches.Patch(color=color, label=label) for label, color in color_regions]

# Add the legend to the plot (you can customize location and other properties)
axs[1].legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.4),
              ncol=4, fontsize=40, frameon=False, title="Genome Regions",title_fontsize=40)
#plt.subplots_adjust(hspace=0.0)  # Adjust the space as needed

plt.suptitle("F-gene mutations frequencies (RSV-A, 2024-2025 season)",
             fontsize=40, fontweight='bold', y=1.02)
fig.savefig("../../plots/all_data/20250321_2429695737_regular_monitoring/nexus_20250321_2429695737_annotated_F_gene_rsva_nonsyn.pdf", format="pdf", bbox_inches="tight")

