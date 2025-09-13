import re
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
import json
import matplotlib.patches as mpatches


with open('../shared/RSV_data_analysis/rsv_definitions/RSVB_nucleotide_mutations_0.9.json', 'r') as file:
    clades_definitions = json.load(file)

european_clades = ["B.D", "B.D.1", "B.D.1.1", "B.D.4", "B.D.4.1", "B.D.4.1.1", "B.D.E.1", "B.D.E.2", "B.D.E.4"]

timeline_tsv_mutation = pd.read_csv(
    '../../preprint/data/timeline_tsv/RSV_B/timeline_mutation_rsv_b_2022_2023_PREPRINT.tsv',
                                    sep='\t',
                                    usecols=['date', 'location','nucleotideMutationFrequency', 'aminoAcidMutationFrequency'])
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
    timeline_tsv_mutation_sorted[str(date)+"_"+location.split(" ")[1]] = sample_muts_sorted
    timeline_tsv_aminoacid_mutation_sorted[str(date)+"_"+location.split(" ")[1]] = sample_aminoacid_muts_sorted



timeline_tsv_mutation_sorted_df = pd.DataFrame.from_dict(timeline_tsv_mutation_sorted)
timeline_tsv_aminoacid_mutation_sorted_df = pd.DataFrame.from_dict(timeline_tsv_aminoacid_mutation_sorted)
df = timeline_tsv_mutation_sorted_df.transpose()
df_aa = timeline_tsv_aminoacid_mutation_sorted_df.transpose()

# drop deletions and insertions:
for mutation in df.columns:

    if (len(re.findall(r'[A-Za-z]+', mutation)[0]) > 1) or (len(re.findall(r'[A-Za-z]+', mutation)[1]) > 1):
        df = df.drop(columns=mutation)
        df_aa = df_aa.loc[:, ~df_aa.columns.str.contains(f'_{re.escape(mutation)}')]

# drop mutations in non-coding regions and synonymous mutations (to save the place in full-genome plot):
df_aa_full_genome_plot = df_aa
df_full_genome_plot = df

to_drop = []
nt_subst_drop = []
for aa_substitution in df_aa_full_genome_plot.columns:

    aa_letters = re.findall(r'[A-Za-z]+', aa_substitution.split('_')[0])
    if len(aa_letters) == 2 and (aa_letters[0] == aa_letters[1]):
        to_drop.append(aa_substitution)
        nt_subst_drop.append( aa_substitution.split('_')[-1])
    elif len(aa_letters) < 2:
        to_drop.append(aa_substitution)
        nt_subst_drop.append( aa_substitution.split('_')[-1])

df_aa_full_genome_plot = df_aa_full_genome_plot.loc[:, ~df_aa_full_genome_plot.columns.isin(to_drop)]
df_full_genome_plot = df_full_genome_plot.loc[:, ~df_full_genome_plot.columns.isin(nt_subst_drop)]


row_separator = next(i for i, label in enumerate(df_full_genome_plot.index) if "GE" in label)

# dataframe of mutation definitions
mutations_df = pd.DataFrame(0, index=european_clades, columns=df_full_genome_plot.columns)
for clade in mutations_df.index:
    for mut in mutations_df.columns:
        mutations_df.at[clade, mut] = 1 if mut in clades_definitions[clade] else 0


sns.set(rc={'figure.figsize': (30, 15)})
sns.set_style("white")
plt.grid(True, linewidth=0.1, color='gray')

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(40, 20), gridspec_kw={'height_ratios': [8, 2]})

sns.heatmap(df_full_genome_plot, ax=axs[0], yticklabels=df_full_genome_plot.index.to_list(),linecolor="black", linewidths=0.0,
            cmap=sns.color_palette("Blues", as_cmap=True),cbar_kws={"shrink": 0.5, "aspect": 10,"label": "Frequency", "pad":0.01})


# Access the colorbar object
colorbar = axs[0].collections[0].colorbar

# Customize the tick labels and label font size
colorbar.ax.tick_params(labelsize=30)  # adjust tick label font size
colorbar.set_label("Frequency", fontsize=35)  # adjust the colorbar label size

sns.heatmap(mutations_df, ax=axs[1], yticklabels=european_clades, linecolor="black", linewidths=0.0,
            cmap=sns.color_palette("Blues", as_cmap=True), cbar_kws={"shrink": 0.5, "pad": 0.01, "aspect": 10})
# colorbar in the second plot is not necessary
colorbar2 = axs[1].collections[0].colorbar
colorbar2.ax.set_visible(False)


# Iterate over the X-axis ticks and color based on "GE" or "ZH"
for ytick in axs[0].get_yticklabels():
    label_text = ytick.get_text()
    if "GE" in label_text:
        ytick.set_color("green")

    elif "ZH" in label_text:
        ytick.set_color("blue")


axs[0].set_facecolor("#ffe6e6")
axs[0].set_xticks([x + 0.5 for x in range(df_full_genome_plot.shape[1])])
axs[0].set_xticklabels([idx.split('_')[0] for idx in df_aa_full_genome_plot.transpose().index], fontsize=25, rotation=90, ha='center', va='top')
axs[0].tick_params(axis='y', labelsize=30)
# Draw a horizontal line to separate ZH and GE
axs[0].hlines(
    y=row_separator,  # Index position where the line should be drawn
    xmin=0,
    xmax=len(df_full_genome_plot.columns),
    colors="black",
    linewidth=2
)
plt.tight_layout()

axs[1].set_xticks([x + 0.5 for x in range(df_full_genome_plot.shape[1])])
axs[1].set_xticklabels(df_full_genome_plot.transpose().index, fontsize=25, rotation=90, ha='center', va='top')
axs[1].set_facecolor("#ffe6e6")
axs[1].tick_params(axis='y', labelsize=27)

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
for xtick in axs[0].get_xticklabels():
    #xtick.set_fontweight("bold")
    xtick_text = xtick.get_text().split('_')[0]
    uppercase_letters = (re.findall(r'[A-Z]', xtick_text))  # Extract genome position

    if len(uppercase_letters) >= 2:  # Ensure there are at least two letters to compare
       # print(uppercase_letters)
        aa_ref = uppercase_letters[0]  # Reference genome position
        aa_alt = uppercase_letters[1]  # Alternative genome position

        # Check if the letters are different
        if aa_ref != aa_alt:
            #xtick.set_color('red')
            xtick.set_fontweight("bold")
for ytick in axs[0].get_yticklabels():
    ytick.set_fontweight("bold")
for ytick in axs[1].get_yticklabels():
    ytick.set_fontweight("bold")

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

# Add the legend to the plot
axs[1].legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.36),
              ncol=5, fontsize=32, frameon=False, title="Genome Regions", title_fontsize=38)
#plt.subplots_adjust(hspace=0.0)  # Adjust the space as needed
plt.suptitle("RSV-B, 2022-2023 season",
             fontsize=30, fontweight='bold', y=1.02)

fig.savefig("../../preprint/plots/mutation_heatmap/timeline_tsv_mutation_sorted_RSVB_22_23_final_revisions.pdf", format="pdf", bbox_inches="tight")

# F-gene heatmap

df_F_gene = df.loc[:, [col for col in df.columns if 5676 < int(re.findall(r'\d+', col)[0]) < 7400]]
df_aa_F_gene = df_aa.loc[:, [col for col in df_aa.columns if 5676 < int(re.findall(r'\d+', col)[-1]) < 7400]]

row_separator = next(i for i, label in enumerate(df_F_gene.index) if "GE" in label)

mutations_df_F_gene = pd.DataFrame(0, index=european_clades, columns=df_F_gene.columns)
for clade in mutations_df_F_gene.index:
    for mut in mutations_df_F_gene.columns:
        mutations_df_F_gene.at[clade, mut] = 1 if mut in clades_definitions[clade] else 0


sns.set_style("white")
plt.grid(True, linewidth=0.1, color='gray')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(40, 20), gridspec_kw={'height_ratios': [8, 2]})

sns.heatmap(df_F_gene, ax=axs[0], yticklabels=df_F_gene.index.to_list(),linecolor="black", linewidths=0.0,
            cmap=sns.color_palette("Blues", as_cmap=True),cbar_kws={"shrink": 0.5, "aspect": 10,"pad": 0.02,"label": "Frequency"})

# Access the colorbar object
colorbar = axs[0].collections[0].colorbar

# Customize the tick labels and label font size
colorbar.ax.tick_params(labelsize=30)  # Adjust tick label font size
colorbar.set_label("Frequency", fontsize=35)  # Adjust the colorbar label size


sns.heatmap(mutations_df_F_gene, ax=axs[1], yticklabels=european_clades, linecolor="black", linewidths=0.0,
            cmap=sns.color_palette("Blues", as_cmap=True), cbar_kws={"shrink": 0.5, "aspect": 10,"pad": 0.02})
# Iterate over the X-axis ticks and color based on "GE" or "ZH"
for ytick in axs[0].get_yticklabels():
    label_text = ytick.get_text()
    if "GE" in label_text:
        ytick.set_color("green")

    elif "ZH" in label_text:
        ytick.set_color("blue")


axs[0].set_facecolor("#ffe6e6")
axs[0].set_xticks([x + 0.5 for x in range(df_F_gene.shape[1])])
axs[0].set_xticklabels([idx.split('_')[0] for idx in df_aa_F_gene.transpose().index], fontsize=35, rotation=90, ha='center', va='top')
axs[0].tick_params(axis='y', labelsize=30)
# Draw a horizontal line to separate ZH and GE
axs[0].hlines(
    y=row_separator,  # Index position where the line should be drawn
    xmin=0,
    xmax=len(df.columns),
    colors="black",
    linewidth=2
)
plt.tight_layout()

bold_regions = [
    (5676, 7400, "F"),  # F region
]

axs[1].set_xticks([x + 0.5 for x in range(df_F_gene.shape[1])])
axs[1].set_xticklabels(df_F_gene.transpose().index, fontsize=35, rotation=90, ha='center', va='top')
axs[1].set_facecolor("#ffe6e6")
axs[1].tick_params(axis='y', labelsize=30)

colorbar2 = axs[1].collections[0].colorbar
colorbar2.ax.set_visible(False)

color_ranges = [
    (5676, 7400, "darkcyan")
]
for xtick in axs[0].get_xticklabels():
    xtick.set_fontweight("bold")

    aa_ref = (re.findall(r'[A-Z]', xtick.get_text())[0])  # Extract genome position
    aa_alt = (re.findall(r'[A-Z]', xtick.get_text())[1])  # Extract genome position

    if aa_ref != aa_alt:
        xtick.set_color('red')

for ytick in axs[0].get_yticklabels():
    ytick.set_fontweight("bold")
for ytick in axs[1].get_yticklabels():
    ytick.set_fontweight("bold")

# Iterate over the X-axis ticks and assign colors based on ranges
for xtick in axs[1].get_xticklabels():
    pos_in_genome = int(re.findall(r'\d+', xtick.get_text())[0])  # Extract genome position
    for start, end, color in color_ranges:
        if start < pos_in_genome < end:
            xtick.set_color(color)
            break  # Stop checking once a match is found
    for start, end, region in bold_regions:
        if start < pos_in_genome < end:

            xtick.set_fontweight("bold")
            break


# Define the color-coded regions as a list of tuples
color_regions = [
    ("F (5676 - 7400)", "darkcyan")
]

# Create a list of mpatches.Patch objects for the legend
legend_patches = [mpatches.Patch(color=color, label=label) for label, color in color_regions]

# Add the legend to the plot
axs[1].legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.5),
              ncol=4, fontsize=35, frameon=False, title="Genome Regions",title_fontsize=35)
#plt.subplots_adjust(hspace=0.0)  # Adjust the space as needed

plt.suptitle("RSV-B, 2022-2023 season",
             fontsize=40, fontweight='bold', y=1.02)
fig.savefig("../../preprint/plots/mutation_heatmap/F_gene_mutations_RSVB_22_23_final_revisions.pdf", format="pdf", bbox_inches="tight")