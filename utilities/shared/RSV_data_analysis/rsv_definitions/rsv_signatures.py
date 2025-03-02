import json
import csv

#Prepare binary data matrix: lineages x signature mutations

#json_file = "RSVB_nucleotide_mutations_0.9.json"
json_file = "RSVA_nucleotide_mutations_0.9.json"
with open(json_file, "r") as file:
    clades_definitions = json.load(file)



# Collect the set of mutations
mutations = set()
for lineage, lineage_mutations in clades_definitions.items():
    mutations.update(lineage_mutations)
# Prepare the CSV data
# First column of lineages, then mutations sorted
header = ['Lineages'] + sorted(mutations, key=lambda x: int(''.join([char for char in x if char.isdigit()])))
rows = []

# Create rows for each lineage with 1 or 0 indicating presence of mutation
for lineage, lineage_mutations in clades_definitions.items():
    row = [lineage] + [1 if mutation in lineage_mutations else 0 for mutation in header[1:]]
    rows.append(row)


# Specify the file name
#filename = 'rsv_b_signatures.csv'
filename = 'rsv_a_signatures.csv'

# Write to CSV
with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(header)  # Write header
    writer.writerows(rows)   # Write data rows

