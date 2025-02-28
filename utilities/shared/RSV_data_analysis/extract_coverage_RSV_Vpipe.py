import pandas as pd
import glob
import argparse

"""
Extract information from 1-based coverage.tsv.gz files after running V-pipe from each sample and
produce a single csv file with coverage
"""

def extract_cov(coverage_tsv_file, sample_name, reference):
    # read coverage.tsv file
    coverage_file = pd.read_csv(coverage_tsv_file, sep = '\t', usecols = ['ref', 'pos', f'{sample_name}/date'])
    # coverage.tsv files are 1-based
    coverage_file = coverage_file[coverage_file['ref']==reference]
    position = pd.DataFrame(coverage_file['pos'])
    coverage = pd.DataFrame(coverage_file[f'{sample_name}/date']).rename(columns={f'{sample_name}/date': 'coverage'})
    total_coverage = pd.concat([position, coverage], axis=1).set_index('pos')
    total_coverage['sample'] = sample_name
    print(total_coverage.head())
    return total_coverage

def process_multiple_coverage_files(input_dir, reference_genome):

    # get list of coverage files in the input directory
    coverage_files = glob.glob(input_dir, recursive = True)
    rows = []
    for coverage_file in coverage_files:
        # extract the sample name from the directory name
        sample_name = coverage_file.split('/')[-4]
        print(sample_name)

        df = extract_cov(coverage_file, sample_name, reference=reference_genome)
        # append the dataframe to the list of dataframes
        rows.append(df)

    coverage_out = pd.concat(rows, axis=0, ignore_index=False)
    return coverage_out


def main(input_dir_cov, reference_genome):
    collected_subsampled_coverage = process_multiple_coverage_files(input_dir_cov, reference_genome)
    collected_subsampled_coverage.to_csv(f'collected_rsv_coverage_{input_dir_cov.split("/")[-6]}.tsv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract information from coverage.tsv.gz files from each sample and produce a single csv file with coverage')
    parser.add_argument("input_dir_cov",
                        help='input directory containing coverage.tsv.gz files')
    parser.add_argument("reference_genome")
    args = parser.parse_args()
    main(args.input_dir_cov, args.reference_genome)