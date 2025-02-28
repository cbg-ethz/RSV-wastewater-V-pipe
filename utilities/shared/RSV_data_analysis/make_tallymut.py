import re
import pandas as pd
import numpy as np
import json
import glob
import argparse
import pysam

"""
The script outputs tallymut.tsv file for variants deconvolution after Lofreq (VCFs)
- extract_cov and process_multiple_coverage_files functions prepares a df of the total coverage for multiple samples (1-based)
- load_convert, process_multiple_vcfs functions process multiple vcf files and make a dataframe of observed mutations x samples (frequency values)
- create_mut_freq_dict function creates a dictionary {'A123C' : 0.67, "T140C" : 0.01, ...} to stores values of frequency from single sample
- make_timeline_mutations_tsv function from each sample collects information into tsv file of columns ["submissionId", "primerProtocol", "reads", "date", "location", "reference", 'nucleotideMutationFrequency']
- make_tallymut_file function prepares a tallymut.tsv file. We iterate through all samples and all signature mutations (of variants used for deconvolution)

Additionally, we set frequency to 0.0 or Nan, depending on coverage (the threshold is pre-specified at the beginning):
- if signature mutation is not called, the value of frequency stays zero if the coverage is sufficient, otherwise - missing value.

"""

COVERAGE_THRESHOLD = 30

# Extract the coverage information from single coverage.tsv.gz file
def extract_cov(coverage_tsv_file, sample_name, reference):
    # read coverage.tsv file
    coverage_file = pd.read_csv(coverage_tsv_file, sep='\t', usecols=['ref', 'pos', f'{sample_name}/date'])
    # coverage.tsv files are 1-based
    coverage_file = coverage_file[coverage_file['ref'] == reference]
    position = pd.DataFrame(coverage_file['pos'])
    coverage = pd.DataFrame(coverage_file[f'{sample_name}/date']).rename(columns={f'{sample_name}/date': 'coverage'})
    total_coverage = pd.concat([position, coverage], axis=1)#.set_index('pos')
    total_coverage['sample'] = sample_name

    return total_coverage

# Extract the coverage information from multiple coverage.tsv.gz files and concatenate them into single data frame
def process_multiple_coverage_files(coverage_input_dir, reference_genome):
    # get list of coverage files in the input directory
    coverage_files = glob.glob(coverage_input_dir, recursive=True)
    rows = []
    for coverage_file in coverage_files:
        # extract the sample name from the directory name
        sample_name = coverage_file.split('/')[-4]
        # print(sample_name)

        df = extract_cov(coverage_file, sample_name, reference=reference_genome)
        # append the dataframe to the list of dataframes
        rows.append(df)

    coverage_out = pd.concat(rows, axis=0, ignore_index=True)
    return coverage_out


# Extract the called mutations and their frequency values from VCF files
def load_convert(vcf_path, sample_name, reference):
    ''' function to load a vcf and output a Pandas dataframe'''
    # Create a VariantFile object
    vcf_file = pysam.VariantFile(vcf_path)
    # record rows for df
    rows = []

    # iterate over the records (variants) in the VCF file
    for record in vcf_file:
        if record.chrom == reference:
            # iterate over possible mutated positions
            for alt in record.alts:
                # Take allele frequences (AF) from INFO
                af = record.info.get('AF')
                # put results in rows
                row = {
                    'sample': sample_name,
                    'pos': record.pos,
                    'ref': record.ref,
                    'alt': alt,
                    'af': af
                }
                rows.append(row)
    # close the VCF file
    vcf_file.close()
    # convert to dataframe
    df_out = pd.DataFrame(rows)

    return df_out


def process_multiple_vcfs(input_dir, reference):
    # get list of VCF files in the input directory
    vcf_files = glob.glob(input_dir, recursive=True)

    # record rows for df
    rows = []
    # iterate over the VCF files from the list
    for vcf_file in vcf_files:
        # extract the sample name from the directory name
        sample_name = vcf_file.split('/')[-5]
        # load the VCF and convert to dataframe
        df = load_convert(vcf_file, sample_name, reference)
        # append the dataframe to the list of rows
        rows.append(df)

    # concatenate all dataframes
    df_out = pd.concat(rows, axis=0, ignore_index=True)

    return (df_out)


def create_mut_freq_dict(col):
    # create empty dictionary to hold the MUTATION PROPORTION FOR SINGLE SAMPLE
    mutation_proportions = {}
    # iterate through all columns except the 'sample' column
    for mutation, frequency in col.items():
        if pd.isna(frequency):
            mutation_proportions[mutation] = None
        else:
            mutation_proportions[mutation] = frequency

    # convert dictionary to a json-formatted string with double quotes
    mutation_json = json.dumps(mutation_proportions, ensure_ascii=False)
    # print(mutation_json)
    return mutation_json


# get frequencies of all observed mutations
def make_timeline_mutations_tsv(path_to_vcf, timeline_tsv, reference):
    # Process multiple vcfs and produce data frame

    output_multiple_vcfs = process_multiple_vcfs(path_to_vcf, reference)
    # print(output_multiple_vcfs)
    output_multiple_vcfs['ref_pos'] = output_multiple_vcfs['ref'] + output_multiple_vcfs['pos'].astype(str)
    output_multiple_vcfs['mut'] = output_multiple_vcfs['ref_pos'] + output_multiple_vcfs['alt']
    output_multiple_vcfs['freq_total'] = output_multiple_vcfs['af']
    output_multiple_vcfs = output_multiple_vcfs[['sample', 'mut', 'freq_total']].drop_duplicates()
    # missing values remains as NaN's
    mut_freq = output_multiple_vcfs.pivot_table(values='freq_total',
                                                index='mut',
                                                columns='sample')

    mut_freq.columns = mut_freq.columns.astype(str)

    tsv_samples_locations = pd.read_csv(timeline_tsv, sep='\t',
                                        usecols=["submissionId", "primerProtocol", "reads", "date", "location", "reference"])

    nucleotide_mut_freq = mut_freq.apply(create_mut_freq_dict, axis=0)

    tsv_samples_locations['nucleotideMutationFrequency'] = tsv_samples_locations["submissionId"].map(
        nucleotide_mut_freq)

    return tsv_samples_locations


def find_loc_code(location):
    code = {
        'Lugano (TI)': '05',
        'Zürich (ZH)': '10',
        'Basel (BS)': '15',
        'Genève (GE)': '16',
        'Chur (GR)': '17',
        'Laupen (BE)': '25'
    }

    return code.get(location, 'Unknown')


def make_tallymut_file(signatures_matrix, timeline_tsv_mutation, collected_coverage):
    rsv_definitions = pd.read_csv(signatures_matrix)
    rsv_definitions.set_index('Lineages', inplace=True)
    signature_muts = rsv_definitions.columns.values


    lineages = rsv_definitions.index
    # columns of tallymut.tsv file
    columns = ['sample',
               'batch',
               'reads',
               'proto',
               'location_code',
               'date',
               'location',
               'pos',
               'gene',
               'base',
               'cov',
               'var',
               'frac'] + list(lineages)

    tallymut = pd.DataFrame(columns=columns)


    mutation_list = []

    # iterate through all SAMPLES:
    for _, row in timeline_tsv_mutation.iterrows():
        nucleotide_mut_freq_raw = row['nucleotideMutationFrequency']

        if pd.isna(nucleotide_mut_freq_raw) or nucleotide_mut_freq_raw == '' or nucleotide_mut_freq_raw is None:
            nucleotide_mut_freq = {}

        else:
            nucleotide_mut_freq = json.loads(row['nucleotideMutationFrequency'])

        location_info = row['location']

        # if sample is complete drop-out -> skip
        if row['submissionId'] not in (collected_coverage['sample'].values):
            print(row['submissionId'])
            continue

        # iterate through all SIGNATURE MUTATIONS:

        for sign_mut in signature_muts:


            sign_pos = int(re.findall(r'\d+', sign_mut)[0])  # position in the genome
            sign_base = re.findall(r'\D+', sign_mut)[1]  # second position = variant nt

            # coverage of the particular sample at the particular position of the signature mutation:
            sign_coverage = (collected_coverage[(collected_coverage['sample'] == row['submissionId']) &
                                                (collected_coverage['pos'] == sign_pos)]['coverage']).iloc[0]

            # if there are no mutations observed in the sample -> set signature frequency to 0.0 or Nan, depending on coverage:
            if sign_coverage > COVERAGE_THRESHOLD:
                sign_freq = 0.0
            else:
                sign_freq = np.nan

            # if mutation is among the keys (observed in the time period of interest), and it's value is not missing (observed in the sample), take the outputted frequency value from Lofreq:

            if (sign_mut in nucleotide_mut_freq.keys()) and (nucleotide_mut_freq[sign_mut] is not None) and sign_freq is not None:
                sign_freq = nucleotide_mut_freq[sign_mut]

            mut_row = {'sample': row['submissionId'],
                       'batch': np.nan,
                       'reads': row['reads'],
                       'proto': row['primerProtocol'],
                       'location_code': find_loc_code(location=location_info),
                       'date': row['date'],
                       'location': location_info,
                       'pos': sign_pos,
                       'gene': np.nan,
                       'base': sign_base,
                       'cov': sign_coverage,
                       'var': np.nan,
                       'frac': sign_freq
                       }

            for lineage in lineages:

                if (rsv_definitions.at[lineage, sign_mut]) == 1:
                    mut_row[lineage] = 'mut'
                else:
                    mut_row[lineage] = np.nan

            mutation_list.append(mut_row)


    # concatenate tallymut.tsv's from separate mutations
    tallymut = pd.concat([tallymut, pd.DataFrame(mutation_list)], ignore_index=True)


    return tallymut


def main(signatures_matrix, input_dir_vcf, input_dir_cov, timeline_tsv, reference_genome):

    # extract mutation frequency information from VCFs and prepare timeline_tsv file - input of make_tallymut_file function
    timeline_tsv_mutation = make_timeline_mutations_tsv(path_to_vcf=input_dir_vcf,
                                                        timeline_tsv=timeline_tsv,
                                                        reference=reference_genome)

    # prepare concat coverage files - input of make_tallymut_file function
    collected_coverage = process_multiple_coverage_files(coverage_input_dir=input_dir_cov,
                                                         reference_genome=reference_genome)

    tallymut = make_tallymut_file(signatures_matrix=signatures_matrix,
                       timeline_tsv_mutation=timeline_tsv_mutation,
                       collected_coverage=collected_coverage)

    tallymut.to_csv(f"tallymut_{reference_genome}.tsv", sep='\t')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='produce tallymut.tsv file, which is used as an input for deconvolution')

    parser.add_argument("signatures_matrix", help='lineages definition matrix with lineages in the rows and mutations in the columns')
    parser.add_argument("input_dir_vcf")
    parser.add_argument("input_dir_cov", help='path to coverage.tsv.gz')
    parser.add_argument("timeline_tsv")
    parser.add_argument("reference_genome")

    args = parser.parse_args()
    main(args.signatures_matrix, args.input_dir_vcf, args.input_dir_cov, args.timeline_tsv, args.reference_genome)
