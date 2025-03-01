import pandas as pd
import requests
import json


url_rsva ='https://lapis.genspectrum.org/rsv-a/sample/aminoAcidMutations?minProportion=0.0&orderBy=sequenceName&limit=10000&downloadAsFile=false'
response_rsva = requests.get(url_rsva)
data_rsva = (json.loads(response_rsva.content))['data']

print("RSV-A mutations")
# Amino acid substitutions of RSV-A observed in the wastewater
rsva_aa_mut = ["F:T12I",
             "F:L20I",
             "F:L20F",
             "F:A23T",
             "F:T91R",
             "F:A103T",
             "F:M115T",
             "F:T122A",
             "F:S377N",
             "F:A518V"]

mutation_proportions_rsva = []
for entry in data_rsva:
    #print(entry)
    if entry['mutation'] in rsva_aa_mut:
        mutation_proportions_rsva.append(
            {
                'mutation':entry['mutation'],
                'percentage %': round(entry['proportion'],5) * 100, # provide in percentage
                'count': entry['count'],
                'total': entry['count'] / entry['proportion']

            }
        )

mutation_proportions_rsva = pd.DataFrame(mutation_proportions_rsva)
print(mutation_proportions_rsva)
#mutation_proportions_rsva.to_csv('rsv_a_globally_found_mutations.csv', index=False)
#print(mutation_proportions)

url_rsvb ='https://lapis.genspectrum.org/rsv-b/sample/aminoAcidMutations?minProportion=0.0&orderBy=sequenceName&limit=10000&downloadAsFile=false'
response_rsvb = requests.get(url_rsvb)
print("RSV-B mutations")
data_rsvb = (json.loads(response_rsvb.content))['data']
#print(data_rsvb)
# Amino acid substitutions of RSV-B observed in the wastewater
rsvb_aa_mut = ["F:F12I",
               "F:R42K",
               "F:S190N",
               "F:R209Q",
               "F:S211N",
               "F:S389P",
               "F:I431L",
               "F:T555A"]



mutation_proportions_rsvb = []
for entry in data_rsvb:
    #print(entry)
    if entry['mutation'] in rsvb_aa_mut:
        mutation_proportions_rsvb.append(
            {
                'mutation':entry['mutation'],
                'percentage % ': round(entry['proportion'],5) * 100,
                'count': entry['count'],
                'total': entry['count']/entry['proportion']
            }
        )

mutation_proportions_rsvb = pd.DataFrame(mutation_proportions_rsvb)
#mutation_proportions_rsvb.to_csv('rsv_b_globally_found_mutations.csv', index=False)

print(mutation_proportions_rsvb)


url_swissb = 'https://lapis.genspectrum.org/rsv-b/sample/aminoAcidMutations?country=Switzerland&minProportion=0.0&orderBy=sequenceName&limit=10000&downloadAsFile=false'
response_rsv_swiss_b = requests.get(url_swissb)
print("swiss RSV-B mutations")
data_rsv_swiss_b = (json.loads(response_rsv_swiss_b.content))['data']


mutation_proportions_rsv_swiss_b = []
for entry in data_rsv_swiss_b:
    #print(entry)
    if entry['mutation'] in rsvb_aa_mut:
        mutation_proportions_rsv_swiss_b.append(
            {
                'mutation':entry['mutation'],
                'percentage % ': round(entry['proportion'],5) * 100,
                'count':entry['count'],
                'total': entry['count']/entry['proportion']
            }
        )

mutation_proportions_rsv_swiss_b = pd.DataFrame(mutation_proportions_rsv_swiss_b)
#mutation_proportions_rsv_swiss_b.to_csv('rsv_b_swiss_found_mutations.csv', index=False)
print(mutation_proportions_rsv_swiss_b)