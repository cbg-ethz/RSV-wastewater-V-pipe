import numpy as np
import pysam
import pandas as pd
import json
import glob
import argparse
import re


THRESHOLD_VALUE = 0.02
RARE_MUTATION_LIMIT_DAYS = 2
COVERAGE_THRESHOLD = 30  # coverage depth below which mutation is treated as missing value


def extract_codon_position(value):
    if pd.isna(value):
        return '*'  # or return None if you'd prefer to have None instead of NaN
    # Extract the first number from the string
    result = int(value)  # re.findall(r'\d+', value)
    if result:
        #    print(result)
        return (result)  # return the first found number as an integer
    return '*'  # if no number is found, return NaN


# Processing vcf files

# Iterate over multiple VCF files and use the name of the directory containing each file as the sample name
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

                CodonPosition = record.info.get('CodonPosition')
                RefAminoAcid = record.info.get('RefAminoAcid')
                AltAminoAcid = record.info.get('AltAminoAcid')
                Gene = record.info.get('Gene')
                # put results in rows
                row = {
                    'sample': sample_name,
                    'pos': record.pos,
                    'ref': record.ref,
                    'alt': alt,
                    'af': af,
                    # 'CodonPosition': extract_codon_position(CodonPosition),
                    'CodonPosition': CodonPosition,
                    'RefAminoAcid': RefAminoAcid[0],
                    'AltAminoAcid': AltAminoAcid[0],
                    'Gene': Gene[0]

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
    # create empty dictionary to hold the mutation proportions for this sample
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


def extract_cov(coverage_tsv_file, sample_name, reference):
    # read coverage.tsv file
    coverage_file = pd.read_csv(coverage_tsv_file, sep='\t', usecols=['ref', 'pos', f'{sample_name}/date'])
    # coverage.tsv files are 1-based
    coverage_file = coverage_file[coverage_file['ref'] == reference]
    position = pd.DataFrame(coverage_file['pos'])
    coverage = pd.DataFrame(coverage_file[f'{sample_name}/date']).rename(columns={f'{sample_name}/date': 'coverage'})
    total_coverage = pd.concat([position, coverage], axis=1).set_index('pos')
    total_coverage['sample'] = sample_name
    # print(total_coverage.head())
    return total_coverage


def process_multiple_coverage_files(input_dir, reference_genome):
    # get list of coverage files in the input directory
    coverage_files = glob.glob(input_dir, recursive=True)
    rows = []
    for coverage_file in coverage_files:
        # extract the sample name from the directory name
        sample_name = coverage_file.split('/')[-4]
        # print(sample_name)

        df = extract_cov(coverage_file, sample_name, reference=reference_genome)
        # append the dataframe to the list of dataframes
        rows.append(df)

    coverage_out = pd.concat(rows, axis=0, ignore_index=False)
    return coverage_out


def main(path_to_vcf, timeline_tsv, path_to_coverage, reference):
    # Process multiple vcfs and produce data frame

    output_multiple_vcfs = process_multiple_vcfs(path_to_vcf, reference)
    # print(output_multiple_vcfs)
    output_multiple_vcfs['ref_pos'] = output_multiple_vcfs['ref'] + output_multiple_vcfs['pos'].astype(str)
    output_multiple_vcfs['mut'] = output_multiple_vcfs['ref_pos'] + output_multiple_vcfs['alt']
    output_multiple_vcfs['freq_total'] = output_multiple_vcfs['af']
    output_multiple_vcfs['AA_mut'] = (output_multiple_vcfs['RefAminoAcid'] +
                                      (output_multiple_vcfs['CodonPosition'].apply(extract_codon_position)).astype(
                                          str) +
                                      output_multiple_vcfs['AltAminoAcid'] +
                                      '_' +
                                      (output_multiple_vcfs['Gene']).astype(str))

    output_multiple_vcfs = output_multiple_vcfs[['sample', 'mut', 'freq_total', 'AA_mut']].drop_duplicates()
    output_multiple_vcfs['mut_aa_mut'] = output_multiple_vcfs['AA_mut'] + '_' + output_multiple_vcfs['mut']
    print(output_multiple_vcfs)
    #    output_multiple_vcfs = output_multiple_vcfs[['sample', 'mut', 'freq_total']].drop_duplicates()
    mut_freq = output_multiple_vcfs.pivot_table(values='freq_total',
                                                index='mut',
                                                columns='sample',
                                                sort=False)

    # For positions where coverage is > threshold, set missing values to zeros (mutation is not present)
    collected_coverage = process_multiple_coverage_files(path_to_coverage, reference)
    collected_coverage = collected_coverage.pivot_table(values='coverage',
                                                        index='pos',
                                                        columns='sample',
                                                        sort=False)
    print(len(mut_freq.index))
    mut_freq.columns = mut_freq.columns.astype(str)

    for mut in mut_freq.index:  # iterate over the mutations
        position = int(re.findall(r'\d+', mut)[0])

        for sample in mut_freq.columns:
            # coverage_at_position = collected_coverage[(collected_coverage.index == position) & (collected_coverage["sample"] == sample)]
            if (collected_coverage.loc[position, sample] >= COVERAGE_THRESHOLD):
                if (pd.isna(mut_freq.loc[mut, sample])):
                    mut_freq.loc[mut, sample] = 0.0
            # if coverage at the position is below the COVERAGE_THRESHOLD -> missing value
            else:
                mut_freq.loc[mut, sample] = np.nan

        # We trust only mutations that appear above THRESHOLD_VALUE for at least RARE_MUTATION_LIMIT_DAYS days
    nonzero_counts = ((mut_freq > THRESHOLD_VALUE) & pd.notna(mut_freq)).sum(axis=1)
    # print(nonzero_counts)
    rare_mutations = nonzero_counts[nonzero_counts < RARE_MUTATION_LIMIT_DAYS]
    # print(rare_mutations)
    print(mut_freq.index.isin(rare_mutations.index))
    # Drop the rows with rare mutations from the DataFrame
    # rare_mutations.to_frame().reset_index(inplace=True)

    #    rare_mutations.index = rare_mutations.index.astype(str)
    # print(rare_mutations)

    #    print(mut_freq[mut_freq.index.isin(rare_mutations.index)])
    #    mut_freq.reset_index(inplace=True)
    print(mut_freq)
    mut_freq_old = mut_freq
    mut_freq = mut_freq[~mut_freq_old.index.isin(rare_mutations.index)]  # ~mut_freq.index.isin(rare_mutations.index)]
    print(len(mut_freq_old.index))
    print(len(rare_mutations.index))
    # set index back
    #    mut_freq.set_index(mut_freq.columns[0])
    #    output_multiple_vcfs = output_multiple_vcfs[~output_multiple_vcfs['mut'].isin(rare_mutations.index)]
    # make missing a.a. mutations unique (otherwise cannot set as an index)

    AA_mut_freq = output_multiple_vcfs.pivot_table(values='freq_total',
                                                   index='mut_aa_mut',
                                                   columns='sample',
                                                   sort=False)

    #    print(keep_mut_numeric_index)
    print(len(AA_mut_freq.index))

    #    AA_mut_freq.reset_index(inplace=True)

    AA_mut_freq = AA_mut_freq.loc[
        ~mut_freq_old.index.isin(rare_mutations.index)]  # ~AA_mut_freq.index.isin(rare_mutations.index)]
    #    AA_mut_freq.set_index(AA_mut_freq.columns[0], inplace=True)

    #    AA_mut_freq = AA_mut_freq[~mut_freq.index.isin(rare_mutations.index)]
    print(AA_mut_freq)
    #    AA_mut_freq.set_index(AA_mut_freq.columns[0], inplace=True)

    #    AA_mut_freq.columns = AA_mut_freq.columns.astype(str) # the order of the columns should be the same as mut_freq

    #    rare_mut_numeric_index = rare_mutations.reset_index().index.tolist()

    tsv_samples_locations = pd.read_csv(timeline_tsv, sep='\t',
                                        usecols=["submissionId", "date", "location", "primerProtocol", "reference"])

    # Each column of mut_freq is passed to the create_mut_freq_dict function as a pandas Series
    # print("new")
    # print(mut_freq.shape)
    # print(AA_mut_freq.shape)
    nucleotide_mut_freq = mut_freq.apply(create_mut_freq_dict, axis=0)
    # print(nucleotide_mut_freq)
    AA_mut_freq_dict = AA_mut_freq.apply(create_mut_freq_dict, axis=0)
    # print(AA_mut_freq_dict)

    tsv_samples_locations['nucleotideMutationFrequency'] = tsv_samples_locations["submissionId"].map(
        nucleotide_mut_freq)

    tsv_samples_locations['aminoAcidMutationFrequency'] = tsv_samples_locations["submissionId"].map(
        AA_mut_freq_dict)

    tsv_samples_locations['lineageFrequencyEstimates'] = None

    tsv_samples_locations.to_csv(
        f'timeline_mutation_{path_to_vcf.split("/")[-7]}.tsv', sep='\t',
        index=False, quoting=3)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process annotated VCF file and tsv files and prepare datamatrix')
    parser.add_argument('path_to_vcf',
                        help='input directory containing annotated VCF files snvs_annotated.vcf')
    parser.add_argument('timeline_tsv', help='path to timeline.tsv file')
    parser.add_argument('path_to_coverage',
                        help='input directory containing coverage tsv files coverage.tsv.gz')
    parser.add_argument('reference',
                        help='reference: it should be e.g. (for RSV-A:) EPI_ISL_412866; (for RSV-B:) EPI_ISL_1653999')

    args = parser.parse_args()

    main(args.path_to_vcf, args.timeline_tsv, args.path_to_coverage, args.reference)