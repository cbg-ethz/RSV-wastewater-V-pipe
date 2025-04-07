#!/usr/bin/env python3

import numpy as np
import pysam
import pandas as pd
import json
import glob
import argparse
import re
from collections import Counter
import yaml
"""
Prepare TSV file in the format which is needed for uploading data to GenSpectrum.
Analyze annotated vcf files.
Only non-synonymous substitutions are included.
"""

THRESHOLD_VALUE = 0.02 # minimal frequency value
RARE_MUTATION_LIMIT_DAYS = 2 # to be included in the heatmap mutation has to appear for at least RARE_MUTATION_LIMIT_DAYS days
COVERAGE_THRESHOLD = 30  # coverage depth below which mutation is treated as missing value


def extract_codon_position(value):
    if pd.isna(value):
        return '*'
    # Extract the first number from the string
    result = int(value)  # re.findall(r'\d+', value)
    return (result)  # return the first found number as an integer

def is_vcf_empty(vcf_path):
    with open(vcf_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                # Found a variant record
                return False
    return True

# Processing vcf files

# Iterate over multiple annotated VCF files and use the name of the directory containing each file as the sample name
def load_convert(vcf_path, sample_name, reference):
    ''' function to load a single  annotated vcf and output a Pandas dataframe'''
    # Create a VariantFile object
    vcf_file = pysam.VariantFile(vcf_path)
    # record rows for df
    rows = []
    for record in vcf_file:
        if record.chrom == reference:
            # iterate over possible mutated positions
            # Get INFO fields defined in the VCF header
            for alt in record.alts:
                # Take allele frequences (AF) from INFO
                af = record.info.get('AF')#
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
                    'Gene': Gene[0]#
                }
                rows.append(row)
    # close the VCF file
    vcf_file.close()
    # convert to dataframe, keys as column names
    df_out = pd.DataFrame(rows)#
    return df_out


def process_multiple_vcfs(input_directories, reference):
    vcf_files_multiple = []
    #print("Entering process_multiple_vcfs function")
    for input_dir in input_directories:
        # get list of VCF files in the input directory
        vcf_files = glob.glob(input_dir, recursive=True)
        vcf_files_multiple.extend(vcf_files)
    
    # record rows for df
    rows = []
    # iterate over the VCF files from the list
    for vcf_file in vcf_files_multiple:
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
    # create empty dictionary to hold the mutation proportions for a particular sample
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

    coverage_file = pd.read_csv(coverage_tsv_file, sep='\t', usecols=[0,1,2])#['ref', 'pos', f'{sample_name}/date'])
    # coverage.tsv files are 1-based
    coverage_file = coverage_file[coverage_file['ref'] == reference]
    position = pd.DataFrame(coverage_file['pos'])
    coverage = pd.DataFrame(coverage_file.iloc[:, 2])
    coverage.columns = ['coverage']
    total_coverage = pd.concat([position, coverage], axis=1).set_index('pos')
    total_coverage['sample'] = sample_name
    # print(total_coverage.head())
    return total_coverage


def process_multiple_coverage_files(input_directories, reference_genome):
    coverage_files_multiple = []
    for input_dir in input_directories:
        # get list of coverage files in the input directory
        coverage_files = glob.glob(input_dir, recursive=True)
        coverage_files_multiple.extend(coverage_files)

    rows = []
    for coverage_file in coverage_files_multiple:
        # extract the sample name from the directory name
        sample_name = coverage_file.split('/')[-4]
        

        df = extract_cov(coverage_file, sample_name, reference=reference_genome)
        # append the dataframe to the list of dataframes
        rows.append(df)

    coverage_out = pd.concat(rows, axis=0, ignore_index=False)
    return coverage_out


def main(path_to_vcf, timeline_tsv, path_to_coverage, reference,path_to_output,BATCH,VIR_STRING):
    if reference == "EPI_ISL_412866":
        subtype = "RSV-A"
    elif reference == "EPI_ISL_1653999":
        subtype = "RSV-B"

    # Process multiple vcfs and produce data frame

    output_multiple_vcfs = process_multiple_vcfs(path_to_vcf, reference)
    # print(output_multiple_vcfs)
    output_multiple_vcfs['ref_pos'] = output_multiple_vcfs['ref'] + output_multiple_vcfs['pos'].astype(str)
    output_multiple_vcfs['mut'] = output_multiple_vcfs['ref_pos'] + output_multiple_vcfs['alt']
    output_multiple_vcfs['freq_total'] = output_multiple_vcfs['af']
    output_multiple_vcfs['AA_mut'] = ((output_multiple_vcfs['Gene']).astype(str) +
                                      ':' +
                                      (output_multiple_vcfs['RefAminoAcid'] +
                                      (output_multiple_vcfs['CodonPosition'].apply(extract_codon_position)).astype(
                                          str) +
                                      output_multiple_vcfs['AltAminoAcid']))


    output_multiple_vcfs = output_multiple_vcfs[['sample', 'mut', 'freq_total', 'AA_mut']].drop_duplicates()

    # Drop deletions and insertions

    to_drop = []
    for nt_substitution in output_multiple_vcfs['mut']:
        if ((len(re.findall(r'[A-Za-z]+', nt_substitution)[0]) > 1) or
                (len(re.findall(r'[A-Za-z]+', nt_substitution)[1]) > 1)):
            to_drop.append(nt_substitution)
    output_multiple_vcfs = output_multiple_vcfs[~output_multiple_vcfs['mut'].isin(to_drop)]

    ###################
    # drop mutations in non-coding regions and synonymous mutations:
    to_drop = []
    for aa_substitution in output_multiple_vcfs['AA_mut']:

        aa_letters = re.findall(r'[A-Za-z]+', aa_substitution.split(':')[1])
        if len(aa_letters) == 2 and (aa_letters[0] == aa_letters[1]):
            to_drop.append(aa_substitution)
        elif len(aa_letters) < 2:
            to_drop.append(aa_substitution)
    output_multiple_vcfs = output_multiple_vcfs[~output_multiple_vcfs['AA_mut'].isin(to_drop)]

    output_multiple_vcfs['mut_aa_mut'] = output_multiple_vcfs['AA_mut']
    output_multiple_vcfs['mut_aa_mut'] = output_multiple_vcfs['AA_mut'] + '_' + output_multiple_vcfs['mut']


    # For positions where coverage is > threshold, set missing values to zeros (mutation is not present,
    # but coverage per position is sufficient)
    collected_coverage = process_multiple_coverage_files(path_to_coverage, reference)
    collected_coverage = collected_coverage.pivot_table(values='coverage',
                                                        index='pos',
                                                        columns='sample',
                                                        sort=False)


    AA_mut_freq = output_multiple_vcfs.pivot_table(values='freq_total',
                                                   index=['AA_mut', 'mut'],
                                                   columns='sample',
                                                   sort=False)

    AA_mut_freq.reset_index(inplace=True)
    no_dupl_aa = []
    #print(AA_mut_freq.head())
    aa_counts = Counter(AA_mut_freq["AA_mut"])

    seen = {}
    for aa in AA_mut_freq["AA_mut"]:
        if aa_counts[aa] > 1:
            #print(aa)
            #print(AA_mut_freq.loc[AA_mut_freq["AA_mut"]==aa])
            count = seen.get(aa, 0) +1
            new_aa = aa + f"_{count}"
            #print(new_aa)
            seen[aa] = count
        else:
            new_aa = aa
        no_dupl_aa.append(new_aa)
    AA_mut_freq["AA_mut"] = no_dupl_aa

    #print(AA_mut_freq.shape)

    mut_freq = AA_mut_freq.drop(columns="AA_mut")
    mut_freq.set_index("mut", inplace=True)
    mut_freq.columns = mut_freq.columns.astype(str)

    for mut in mut_freq.index:  # iterate over the mutations
        position = int(re.findall(r'\d+', mut)[0])

        for sample in mut_freq.columns:

            if (collected_coverage.loc[position, sample] >= COVERAGE_THRESHOLD):
                if pd.isna(mut_freq.loc[mut, sample]):
                    mut_freq.loc[mut, sample] = 0.0
            # if coverage at the position is below the COVERAGE_THRESHOLD -> missing value
            else:
                mut_freq.loc[mut, sample] = np.nan

    # We keep only mutations that appear above THRESHOLD_VALUE for at least RARE_MUTATION_LIMIT_DAYS days
    nonzero_counts = ((mut_freq > THRESHOLD_VALUE) & pd.notna(mut_freq)).sum(axis=1)

    rare_mutations = nonzero_counts[nonzero_counts < RARE_MUTATION_LIMIT_DAYS]


    AA_mut_freq = AA_mut_freq.loc[
        ~AA_mut_freq["mut"].isin(rare_mutations.index)]

    AA_mut_freq_sorted = AA_mut_freq.sort_values(
        by="mut",
        key=lambda x: x.apply(lambda y: int(re.findall(r'\d+', y)[0]))
    )
    AA_mut_freq = AA_mut_freq_sorted

    for mut in AA_mut_freq["mut"]:  # iterate over the mutations
        position = int(re.findall(r'\d+', mut)[0])
        for sample in AA_mut_freq.columns[2:]: # skip columns 'mut' and 'AA_mut'
            index_AA = AA_mut_freq.index[AA_mut_freq["mut"]==mut][0]
            if (collected_coverage.at[position, sample] >= COVERAGE_THRESHOLD):
                if pd.isna(AA_mut_freq.loc[index_AA, sample]):
                    AA_mut_freq.loc[index_AA, sample] = 0.0
            # if coverage at the position is below the COVERAGE_THRESHOLD -> missing value
            else:
                AA_mut_freq.loc[index_AA, sample] = np.nan

    AA_dataframe = AA_mut_freq.drop(columns="mut")
    AA_dataframe.set_index("AA_mut", inplace=True)
    nt_dataframe = AA_mut_freq.drop(columns="AA_mut")
    nt_dataframe.set_index("mut", inplace=True)

    tsv_samples_locations = pd.read_csv(timeline_tsv, sep='\t',
                                        usecols=["submissionId", "date", "location", "primerProtocol", "reference"])

    nucleotide_mut_freq = nt_dataframe.apply(create_mut_freq_dict, axis=0)
    # print(nucleotide_mut_freq)
    AA_mut_freq_dict = AA_dataframe.apply(create_mut_freq_dict, axis=0)
    # print(AA_mut_freq_dict)

    tsv_samples_locations['nucleotideMutationFrequency'] = tsv_samples_locations["submissionId"].map(
        nucleotide_mut_freq)

    tsv_samples_locations['aminoAcidMutationFrequency'] = tsv_samples_locations["submissionId"].map(
        AA_mut_freq_dict)

    tsv_samples_locations['lineageFrequencyEstimates'] = None
######### Add the search for possible substrings ##########
    # Remove the subtype that is not the reference type
    #tsv_samples_locations = tsv_samples_locations[(tsv_samples_locations['reference']== subtype)
    #                                              | (tsv_samples_locations['reference']=='RSV A and B')
    #                                              | (tsv_samples_locations['reference']=='RSV subtype A and B')
    #                                              | (tsv_samples_locations['reference']=='sequencing of both RSV-A and B')
    #                                              | (tsv_samples_locations['reference']=='RSV-A and B')]
    
    unique_refs=np.unique(tsv_samples_locations['reference'])
    #print(unique_refs)
    possible_vir_strings=VIR_STRING.split("|")
    #print(possible_vir_strings)
    # Find items in first_list that are NOT in second_list
    not_in_second = [s for s in unique_refs if s not in possible_vir_strings]
    if len(not_in_second) > 0:
        print("WARNING: the following references is not is the provided refstring:", not_in_second)

    tsv_samples_locations = tsv_samples_locations[tsv_samples_locations['reference'].apply(lambda x: (subtype in possible_vir_strings))]
    ####

    tsv_samples_locations["reference"] = subtype

    #tsv_samples_locations.to_csv(
    #    f'timeline_mutation_{reference}_annotated_nonsyn.tsv', sep='\t',
    #    index=False, quoting=3)
###### adjusted output format #######
    tsv_samples_locations.to_csv(
    f'{path_to_output}{BATCH}_{reference}_Mutations_Dashboard.tsv', sep='\t',
    index=False, quoting=3)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process VCF file and tsv files and prepare datamatrix')
    parser.add_argument('--vpipe_dir',  # added vpipe base dir as input 
                        help='basedir to vpipe. e.g. /cluster/project/pangolin/rsv_pipeline/working')
    parser.add_argument('--path_to_vcf', nargs='+', 
                        help='input directory containing VCF files')
    parser.add_argument('--timeline_tsv', help='path to timeline.tsv file')
    parser.add_argument('--path_to_coverage', nargs='+',
                        help='input directory containing coverage tsv files')
    parser.add_argument('--reference',
                        help='reference: it should be e.g. (for RSV-A:) EPI_ISL_412866; (for RSV-B:) EPI_ISL_1653999')
    parser.add_argument('--config',
                        help='path to config file to read from') #added input of config file to read from
    parser.add_argument('--virus_string',
                        help='all possible accepted string for the specific virus') #added input of config file to read from
    
    args = parser.parse_args()

    # Read the YAML file
    with open(args.config, 'r') as file:
        configs = yaml.safe_load(file)  # Read and parse the YAML file  
    
    # in the samples .tsv file the BATCH is saved
    # we want to get the batch name so that we can use it in the output file name
    samples_tsv = configs['input']['samples_file']
   
    # Read the TSV file into a DataFrame
    samples_tsv_df = pd.read_csv(samples_tsv, sep="\t")   

    # Ensure the second column exists and calculate the max value
    if len(samples_tsv_df.columns) > 1:
        latest_batch = samples_tsv_df.iloc[:, 1].max()  # Access the second column using iloc, latest batch is accessed with max()
    else:
        print("The samples.tsv file does not have a second column to extract the latest batch.")

    out_dir=f'{args.vpipe_dir}/MutationFrequencies/'

    main(args.path_to_vcf, args.timeline_tsv, args.path_to_coverage, args.reference, out_dir, latest_batch, args.virus_string)
