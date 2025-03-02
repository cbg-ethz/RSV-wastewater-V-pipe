import requests
import json

# To call mutation lineage-defining we set threshold of frequency at 0.9 (default).
proportion = 0.90
# """
# Extract RSV-B clades definitions:
# """
clades = [

    "B.D",
    "B.D.1",
    "B.D.1.1",
    "B.D.2",
    "B.D.3",
    "B.D.4",
    "B.D.4.1",
    "B.D.4.1.1",
    "B.D.E.1",
    "B.D.E.2",
    "B.D.E.3",
    "B.D.E.4"
]

# RSV-B nucleotide mutations
clades_definitions = {}
# Create a list of mutations for each clade:
for clade in clades:
    url = f'https://lapis.genspectrum.org/rsv-b/sample/nucleotideMutations?lineage={clade}&minProportion={proportion}&limit=1000&dataFormat=JSON&downloadAsFile=false'
    response = requests.get(url)
    data = (json.loads(response.content))
   # print(data)
    mutations = [single_mut_info["mutation"] for single_mut_info in data['data']]
    clades_definitions[clade] = mutations

print(clades_definitions)
file_path = f"RSVB_nucleotide_mutations_{proportion}.json"
with open(file_path, 'w') as file:
    json.dump(clades_definitions, file, indent=4)


# """
# Extract RSV-A clades definitions:
# """

clades = [
"A.D",

"A.D.1",
"A.D.1.1",
"A.D.1.2",
"A.D.1.3",
"A.D.1.4",
"A.D.1.5",
"A.D.1.6",
"A.D.1.7",
"A.D.1.8",

"A.D.2",
"A.D.2.1",
"A.D.2.2",
"A.D.2.2.1",

"A.D.3",
"A.D.3.1",
"A.D.3.2",
"A.D.3.3",
"A.D.3.4",
"A.D.3.5",
"A.D.3.6",

"A.D.4",
"A.D.4.1",

"A.D.5",
"A.D.5.1",
"A.D.5.2",
"A.D.5.3",
"A.D.5.4"
]

# RSV-A nucleotide mutations
clades_definitions = {}
# Create a list of mutations for each clade:
for clade in clades:
    url = f'https://lapis.genspectrum.org/rsv-a/sample/nucleotideMutations?lineage={clade}&minProportion={proportion}&limit=1000&dataFormat=JSON&downloadAsFile=false'
    response = requests.get(url)
    data = (json.loads(response.content))
   # print(data)
    mutations = [single_mut_info["mutation"] for single_mut_info in data['data']]
    clades_definitions[clade] = mutations

print(clades_definitions)
file_path = f"RSVA_nucleotide_mutations_{proportion}.json"
with open(file_path, 'w') as file:
    json.dump(clades_definitions, file, indent=4)

