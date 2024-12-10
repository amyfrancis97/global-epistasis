import numpy as np
import pandas as pd
import os
import argparse
import mavenn

def get_absolut_input(absolut_final_bindings_path, mutant_file_path, output_dir, antigen, cdr_selection, percentage_prefix, filter_epitope=False):
    """
    Converts Absolut! output into a file compatible with MAVE-NN input.

    Parameters:
    absolut_final_bindings_path (str): The file and location of the Absolut! output of bindings.
    mutant_file_path (str): The file and location of the mutant file containing the CDR mutants and the hamming distances generated from the '1_get_mutants.py' script.
    output_dir (str): Specified location in which to save the output files.
    antigen (str): Antigen name to be used as file naming prefix.
    cdr_selection (str): CDR selection for file naming.
    filter_epitope (str or bool): Epitope to filter rows by, or False if no filtering.
    
    Returns:
    None: Outputs MAVE-NN-compatible files for both all_mutants and single_epitope (if applicable).
    """
        
    # Read binding data from absolut
    binding_data = pd.read_csv(absolut_final_bindings_path, sep="\t", skiprows=1)
    binding_data = binding_data.rename(columns={'CDR3': 'x'})

    # Read in original CDR sequences with number of mutants per sequence
    CDR_data = pd.read_csv(mutant_file_path, sep="\t")

    # Merge the datasets
    all_mutants_restructured = CDR_data.merge(binding_data, how='outer')[['dist', 'Energy', 'x', 'Structure']].dropna().reset_index(drop=True)

    # Create a filtered version if filter_epitope is provided
    if filter_epitope:
        single_epitope_data = all_mutants_restructured[all_mutants_restructured['Structure'] == filter_epitope].drop("Structure", axis=1)
    else:
        single_epitope_data = None  # No filtered data if no epitope is specified

    all_mutants_restructured = all_mutants_restructured.drop("Structure", axis=1)  # Remove Structure column for full dataset

    # Function to split and save data
    def split_and_save(data, prefix):
        # Assign 'training', 'validation', 'test' to non-single mutants (dist != 1)
        mask_non_single = data['dist'] != 1
        data.loc[mask_non_single, 'set'] = np.random.choice(
            a=['training', 'validation', 'test'],
            p=[0.6, 0.2, 0.2],
            size=mask_non_single.sum()
        )

        # Explicitly assign 'training' to all single mutants (dist == 1)
        data.loc[data['dist'] == 1, 'set'] = 'training'
    
        # Replicate y to dy if no sigma available
        data.insert(3, 'dy', data['Energy'])
    
        data = data.rename(columns={'dist': 'dist', 'Energy': 'y'})
        data = data.drop_duplicates(subset="x", keep="first")

        # Split the dataset for MAVE-NN
        trainval_df, test_df = mavenn.split_dataset(data)

        # Save the trainval and test datasets
        trainval_df.to_csv(f'{output_dir}/{prefix}_{percentage_prefix}_trainval_df.txt', sep="\t", index=None)
        test_df.to_csv(f'{output_dir}/{prefix}_{percentage_prefix}_test_df.txt', sep="\t", index=None)

    # Save all mutants data
    split_and_save(all_mutants_restructured, f'{antigen}_{cdr_selection}_all_mutants')

    # Save single epitope data if filter_epitope is specified
    if single_epitope_data is not None:
        split_and_save(single_epitope_data, f'{antigen}_{cdr_selection}_single_epitope')

# %%
def main():
    parser = argparse.ArgumentParser(description="Get MAVE-NN input from Absolut! output.")
    parser.add_argument('--absolut_final_bindings_path', type=str, help="Absolute final bindings file location.")
    parser.add_argument('--mutant_file_path_with_dist', type=str, help="Mutant file path containing distances generated from '1_get_mutants.py'.")
    parser.add_argument('--output_dir', type=str, default='.', help="Directory to save output files.")
    parser.add_argument('--antigen', type=str, default='.', help="Name extension to save output files.")
    parser.add_argument('--cdr_selection', type=str, default='.', help="Name extension to save output files.")
    parser.add_argument('--filter_epitope', type=str, default=False, help="Filter rows by epitope, defaults to False (no filtering).")
    parser.add_argument('--percentage_prefix', type=str, default='.', help="Name extension for percentages.")

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Call get_absolut_input with the filter_epitope argument
    get_absolut_input(
        args.absolut_final_bindings_path,
        args.mutant_file_path_with_dist,
        args.output_dir,
        args.antigen,
        args.cdr_selection,
        args.percentage_prefix,
        args.filter_epitope
    )

# %%
if __name__ == '__main__':
    main()
# %%

